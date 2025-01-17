#include <string>
#include <stdexcept>
#include <fstream>
#include <fmt/format.h>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>

#include "chargefw2.h"
#include "mmcif.h"
#include "common.h"
#include "bonds.h"
#include "../config.h"
#include "../structures/bond.h"
#include "../periodic_table.h"
#include "../utility/strings.h"
#include "../exceptions/file_exception.h"


void mmCIF::read_protein_molecule(gemmi::cif::Block &data, std::unique_ptr<std::vector<Atom>> &atoms) {
    auto structure = gemmi::make_structure_from_block(data);
    if (structure.models.empty()) {
        throw std::runtime_error("Not enough information to create a structure");
    }

    /* Read first model only */
    auto model = structure.models[0];
    size_t idx = 0;
    for (const auto &chain: model.chains) {
        for (const auto &residue: chain.residues) {
            bool hetatm = residue.het_flag == 'H';
            for (const auto &atom: residue.atoms) {
                double x = atom.pos.x;
                double y = atom.pos.y;
                double z = atom.pos.z;
                auto element = PeriodicTable::pte().get_element_by_symbol(get_element_symbol(atom.element.name()));

                if(atom.charge) {
                    fmt::print("Got charge {} on{}\n", atom.charge, atom.element.name());
                }

                if (not atom.has_altloc() or atom.altloc == 'A') {
                    if ((not hetatm) or
                        (config::read_hetatm and residue.name != "HOH") or
                        (config::read_hetatm and not config::ignore_water)) {
                        atoms->emplace_back(idx, element, x, y, z, atom.name, residue.seqid.num.value, residue.name, chain.name, hetatm);
                        atoms->back()._set_formal_charge(atom.charge);
                        idx++;
                    }
                }
            }
        }
    }
}


void mmCIF::read_ccd_molecule(gemmi::cif::Block &data, std::unique_ptr<std::vector<Atom>> &atoms, std::unique_ptr<std::vector<Bond>> &bonds) {
    auto atom_table = data.find("_chem_comp_atom.", {"model_Cartn_x", "model_Cartn_y", "model_Cartn_z",
                                                     "type_symbol", "atom_id", "comp_id", "charge"});
    size_t idx = 0;

    std::map<std::string, const Atom*> atom_names;
    for (const auto row: atom_table) {
        double x, y, z;
        try {
            x = std::stod(row[0]);
            y = std::stod(row[1]);
            z = std::stod(row[2]);
        } catch (std::exception &) {
            throw std::runtime_error("Cannot load coordinates");
        }

        auto element = PeriodicTable::pte().get_element_by_symbol(get_element_symbol(row[3]));
        auto atom_name = row[4];
        auto residue = row[5];
        int charge  = 0;
        try {
            charge = std::stoi(row[6]);
        } catch (std::exception &){
            /* Keep default */
        }
        auto residue_id = 0;

        if (charge) {
            fmt::print("Got charge {} on {}\n", charge, atom_name);
        }

        atoms->emplace_back(idx, element, x, y, z, atom_name, residue_id, residue, "", false);
        atoms->back()._set_formal_charge(charge);
        idx++;
    }

    for(const auto &atom: *atoms) {
        atom_names[atom.name()] = &atom;
    }

    auto bond_table = data.find("_chem_comp_bond.", {"atom_id_1", "atom_id_2", "value_order"});
    for (const auto row: bond_table) {
        std::string atom1_name = row[0];
        std::string atom2_ = row[1];
        std::string order_str = row[2];
        int order;
        if (order_str == "SING") {
            order = 1;
        } else if (order_str == "DOUB") {
            order = 2;
        } else if (order_str == "TRIP") {
            order = 3;
        } else {
            continue;
        }
        bonds->emplace_back(atom_names[row[0]], atom_names[row[1]], order);
    }
}


void mmCIF::process_record(const std::string &structure_data, std::unique_ptr<std::vector<Molecule>> &molecules) {

    auto atoms = std::make_unique<std::vector<Atom>>();
    auto bonds = std::make_unique<std::vector<Bond>>();

    std::string name;
    try {
        gemmi::cif::Document doc = gemmi::cif::read_string(structure_data);
        auto data = doc.sole_block();
        name = data.name;
        const auto names = data.get_mmcif_category_names();

        bool has_atom_site = std::find(names.begin(), names.end(), "_atom_site.") != names.end();
        bool has_chem_comp = std::find(names.begin(), names.end(), "_chem_comp_atom.") != names.end();

        if (has_atom_site) {
            read_protein_molecule(data, atoms);
        } else if (has_chem_comp) {
            read_ccd_molecule(data, atoms, bonds);
        }

        if (atoms->empty()) {
            throw std::runtime_error("No atoms were loaded");
        }

        if (has_atom_site) {
            bonds = get_bonds(atoms);
        }

        molecules->emplace_back(name, std::move(atoms), std::move(bonds));

    } catch (std::exception &e) {
        fmt::print(stderr, "Error when reading {}: {}\n", name, e.what());
    }
}


MoleculeSet mmCIF::read_file(const std::string &filename) {

    std::string line;
    std::string structure_data;

    auto molecules = std::make_unique<std::vector<Molecule>>();
    try {
        std::ifstream file(filename);
        if (!file) {
            throw FileException("Cannot open file: " + filename);
        }
        while(std::getline(file, line)) {
            if (line.starts_with("#") or line.empty()) {
                continue;
            }
            if (line.starts_with("data_")) {
                if (not structure_data.empty()) {
                    process_record(structure_data, molecules);
                }
                structure_data = line;
            } else {
                structure_data += "\n" + line;
            }
        }
        if (structure_data.empty()) {
            throw std::runtime_error("Empty record");
        }
        process_record(structure_data, molecules);
    }
    catch (std::exception &e) {
        throw FileException("Cannot load structure from file: " + filename + e.what());
    }
    return MoleculeSet(std::move(molecules));
}

mmCIF::mmCIF() = default;
