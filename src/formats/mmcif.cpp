//
// Created by krab1k on 28.1.19.
//

#include <string>
#include <map>
#include <set>
#include <sstream>
#include <fstream>
#include <fmt/format.h>
#include <boost/algorithm/string.hpp>

#include "chargefw2.h"
#include "../config.h"
#include "mmcif.h"
#include "common.h"
#include "bonds.h"
#include "../periodic_table.h"


const std::set<std::string> required_protein_fields = {"Cartn_x", "Cartn_y", "Cartn_z", "type_symbol",
                                                       "label_atom_id", "label_alt_id", "label_seq_id",
                                                       "auth_seq_id", "label_comp_id", "label_asym_id",
                                                       "pdbx_PDB_model_num"};

const std::set<std::string> required_ccd_fields = {"model_Cartn_x", "model_Cartn_y", "model_Cartn_z", "type_symbol",
                                                   "atom_id", "comp_id"};


bool has_required_keys(const std::map<std::string, size_t> &fields, const std::set<std::string> &required,
                       std::string &missing);

void
read_protein_molecule(std::ifstream &file, const std::string &name, std::unique_ptr<std::vector<Molecule>> &molecules);

void read_ccd_molecule(std::ifstream &file, const std::string &name, std::unique_ptr<std::vector<Molecule>> &molecules);


MoleculeSet mmCIF::read_file(const std::string &filename) {
    std::ifstream file(filename);
    if (!file) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    auto molecules = std::make_unique<std::vector<Molecule> >();

    std::string line;
    std::string name;
    try {
        while (std::getline(file, line)) {
            if (boost::starts_with(line, "_entry.id")) {
                std::stringstream ss(line);
                ss >> name >> name;
                read_protein_molecule(file, name, molecules);
                break;
            }

            if (boost::starts_with(line, "_chem_comp.id")) {
                std::stringstream ss(line);
                ss >> name >> name;
                read_ccd_molecule(file, name, molecules);
                break;
            }
        }
    } catch (const std::exception &) {
        fmt::print(stderr, "Invalid mmCIF file\n");
        exit(EXIT_FILE_ERROR);
    }

    return MoleculeSet(std::move(molecules));
}


bool has_required_keys(const std::map<std::string, size_t> &fields, const std::set<std::string> &required,
                       std::string &missing) {
    for (const auto &req_field: required) {
        if (fields.find(req_field) == fields.end()) {
            missing = req_field;
            return false;
        }
    }
    return true;
}


void
read_protein_molecule(std::ifstream &file, const std::string &name, std::unique_ptr<std::vector<Molecule>> &molecules) {

    std::string line;
    while (std::getline(file, line)) {
        if (boost::starts_with(line, "_atom_site.")) {
            break;
        }
    }

    auto atoms = std::make_unique<std::vector<Atom>>();
    std::map<std::string, size_t> record_positions;
    size_t category_idx = 0;
    do {
        auto it = line.find('.');
        auto category = line.substr(it + 1);
        boost::trim(category);

        record_positions[category] = category_idx;
        category_idx++;
        std::getline(file, line);
    } while (boost::starts_with(line, "_atom_site."));

    std::string missing;
    if (not has_required_keys(record_positions, required_protein_fields, missing)) {
        fmt::print(stderr, "mmCIF file does not contain required field: _atom_site.{}\n", missing);
        exit(EXIT_FILE_ERROR);
    }

    size_t idx = 0;

    int first_model_no = -1;
    do {
        std::stringstream ss(line);
        std::vector<std::string> records{std::istream_iterator<std::string>(ss),
                                         std::istream_iterator<std::string>{}};

        bool hetatm = false;
        if (line[0] == 'H') {
            hetatm = true;
        }

        auto it = record_positions.find("Cartn_x");
        auto x = std::stod(records[it->second]);

        it = record_positions.find("Cartn_y");
        auto y = std::stod(records[it->second]);

        it = record_positions.find("Cartn_z");
        auto z = std::stod(records[it->second]);

        it = record_positions.find("type_symbol");
        auto symbol = records[it->second];
        auto element = PeriodicTable::pte().get_element_by_symbol(get_element_symbol(symbol));

        it = record_positions.find("label_atom_id");
        auto atom_name = fix_atom_name(records[it->second]);

        it = record_positions.find("label_alt_id");
        auto alt_id = records[it->second];

        int residue_id;
        /* HETATM record does not have label_seq_id */
        if (not hetatm) {
            it = record_positions.find("label_seq_id");
            residue_id = std::stoi(records[it->second]);
        } else {
            it = record_positions.find("auth_seq_id");
            residue_id = std::stoi(records[it->second]);
        }

        it = record_positions.find("label_comp_id");
        auto residue = records[it->second];

        it = record_positions.find("label_asym_id");
        auto chain_id = records[it->second];

        it = record_positions.find("pdbx_PDB_model_num");
        auto model_no = std::stoi(records[it->second]);

        if (first_model_no == -1) {
            first_model_no = model_no;
        } else {
            if (first_model_no < model_no) {
                /* We are currently reading atom from another model */
                break;
            }
        }

        if (alt_id == "." or not is_already_loaded(*atoms, atom_name, residue_id)) {
            if ((not hetatm) or
                (config::read_hetatm and residue != "HOH") or
                (config::read_hetatm and not config::ignore_water)) {
                atoms->emplace_back(idx, element, x, y, z, atom_name, residue_id, residue, chain_id, hetatm);
                idx++;
            }
        }

        std::getline(file, line);
        boost::trim(line);
    } while (line[0] != '#');

    auto bonds = get_bonds(atoms);
    std::map<size_t, int> charges;
    molecules->emplace_back(name, std::move(atoms), std::move(bonds), charges);
}


void
read_ccd_molecule(std::ifstream &file, const std::string &name, std::unique_ptr<std::vector<Molecule>> &molecules) {
    std::string line;
    while (std::getline(file, line)) {
        if (boost::starts_with(line, "_chem_comp_atom.")) {
            break;
        }
    }

    auto atoms = std::make_unique<std::vector<Atom>>();
    std::map<std::string, size_t> record_positions;
    size_t category_idx = 0;
    do {
        auto it = line.find('.');
        auto category = line.substr(it + 1);
        boost::trim(category);

        record_positions[category] = category_idx;
        category_idx++;
        std::getline(file, line);
    } while (boost::starts_with(line, "_chem_comp_atom."));

    std::string missing;
    if (not has_required_keys(record_positions, required_ccd_fields, missing)) {
        fmt::print(stderr, "mmCIF file does not contain required field: _chem_comp_atom.{}\n", missing);
        exit(EXIT_FILE_ERROR);
    }

    size_t idx = 0;
    do {
        std::stringstream ss(line);
        std::vector<std::string> records{std::istream_iterator<std::string>(ss),
                                         std::istream_iterator<std::string>{}};

        auto it = record_positions.find("model_Cartn_x");
        auto x = std::stod(records[it->second]);

        it = record_positions.find("model_Cartn_y");
        auto y = std::stod(records[it->second]);

        it = record_positions.find("model_Cartn_z");
        auto z = std::stod(records[it->second]);

        it = record_positions.find("type_symbol");
        auto symbol = records[it->second];
        auto element = PeriodicTable::pte().get_element_by_symbol(get_element_symbol(symbol));

        it = record_positions.find("atom_id");
        auto atom_name = fix_atom_name(records[it->second]);

        int residue_id = 0;

        it = record_positions.find("comp_id");
        auto residue = records[it->second];

        atoms->emplace_back(idx, element, x, y, z, atom_name, residue_id, residue, "0", false);
        idx++;

        std::getline(file, line);
        boost::trim(line);
    } while (line[0] != '#');

    auto bonds = get_bonds(atoms);
    std::map<size_t, int> charges;
    molecules->emplace_back(name, std::move(atoms), std::move(bonds), charges);
}
