#include <string>
#include <stdexcept>
#include <fstream>
#include <format>
#include <gemmi/read_cif.hpp>
#include <gemmi/mmcif.hpp>

#include "mmcif.h"
#include "common.h"
#include "bonds.h"
#include "../structures/bond.h"
#include "../periodic_table.h"
#include "../utility/exceptions.h"


void mmCIF::read_protein_molecule(gemmi::cif::Block &data, std::unique_ptr<std::vector<Atom>> &atoms) {
    auto structure = gemmi::make_structure_from_block(data);
    if (structure.models.empty()) {
        throw std::runtime_error("The provided structure has no models");
    }

    /* Read the first model only */
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

                if (keep_atom(atom, residue)) {
                    atoms->emplace_back(idx, element, x, y, z, atom.name, residue.seqid.num.value, residue.name, chain.name, hetatm);
                    atoms->back()._set_formal_charge(atom.charge);
                    idx++;
                }
            }
        }
    }
}

void mmCIF::process_record(const std::string &structure_data, std::unique_ptr<std::vector<Molecule>> &molecules) {

    auto atoms = std::make_unique<std::vector<Atom>>();
    auto bonds = std::make_unique<std::vector<Bond>>();

    gemmi::cif::Document doc = gemmi::cif::read_string(structure_data);
    auto data = doc.sole_block();
    std::string name = data.name;
    const auto names = data.get_mmcif_category_names();

    if (std::ranges::find(names, "_atom_site.") != names.end()) {
        read_protein_molecule(data, atoms);
    }
    else {
        throw std::runtime_error("The mmCIF file does not have _atom_site category");
    }

    if (atoms->empty()) {
        throw std::runtime_error("No atoms were loaded");
    }

    bonds = get_bonds(atoms);
    molecules->emplace_back(name, std::move(atoms), std::move(bonds));
}


MoleculeSet mmCIF::read_file(const std::string &filename) {

    std::string line;
    std::string structure_data;

    auto molecules = std::make_unique<std::vector<Molecule>>();
    try {
        std::ifstream file(filename);
        if (!file) {
            throw FileException(std::format("Cannot open file: {}", filename));
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
        throw FileException(std::format("Cannot load structure: {}", e.what()));
    }
    return MoleculeSet(std::move(molecules));
}

mmCIF::mmCIF() = default;
