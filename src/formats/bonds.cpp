#include <vector>
#include <memory>
#include <string>
#include <tuple>
#include <fstream>
#include <filesystem>
#include <map>
#include <sstream>

#include "chargefw2.h"
#include "bonds.h"
#include "common.h"
#include "../structures/atom.h"
#include "../structures/bond.h"
#include "../exceptions/internal_exception.h"

namespace fs = std::filesystem;


void load_residues_info(const std::string &filename,
                        std::map<std::string, std::vector<std::tuple<std::string, std::string, int>>> &residues_data);


void update_bonds(std::unique_ptr<std::vector<Bond>> &bonds, const std::map<std::string, const Atom *> &residue_atoms);


void load_residues_info(const std::string &filename,
                        std::map<std::string, std::vector<std::tuple<std::string, std::string, int>>> &residues_data) {
    std::ifstream file(filename);
    if (!file) {
        throw InternalException(fmt::format("Unable to open amino acids data file: {}", filename));
    }

    std::string line;
    while (std::getline(file, line)) {
        std::string residue = line;
        std::getline(file, line);
        while (not line.empty()) {
            std::string atom1_name;
            std::string atom2_name;
            int bond_order = -1;

            std::stringstream ss(line);
            ss >> atom1_name >> atom2_name >> bond_order;

            atom1_name = fix_atom_name(atom1_name);
            atom2_name = fix_atom_name(atom2_name);

            residues_data[residue].emplace_back(atom1_name, atom2_name, bond_order);
            std::getline(file, line);
        }
    }
}


void update_bonds(std::unique_ptr<std::vector<Bond>> &bonds, const std::map<std::string, const Atom *> &residue_atoms) {

    static std::map<std::string, std::vector<std::tuple<std::string, std::string, int>>> residues_data;
    static bool basic_loaded = false;
    static bool all_loaded = false;

    if (not basic_loaded) {
        load_residues_info(fs::path(INSTALL_DIR) / "share" / "amino_acids.txt", residues_data);
        basic_loaded = true;
    }

    auto residue = residue_atoms.begin()->second->residue();
    auto it = residues_data.find(residue);

    /* Try to load additional residues */
    if (not all_loaded and it == residues_data.end()) {
        load_residues_info(fs::path(INSTALL_DIR) / "share" / "other_residues.txt", residues_data);
        all_loaded = true;
        it = residues_data.find(residue);
    }

    if (it != residues_data.end()) {
        for (const auto &[atom1_name, atom2_name, order]: it->second) {
            auto it1 = residue_atoms.find(atom1_name);
            auto it2 = residue_atoms.find(atom2_name);

            if (it1 != residue_atoms.end() and it2 != residue_atoms.end()) {
                bonds->emplace_back(it1->second, it2->second, order);
            }
        }
    }
}


std::unique_ptr<std::vector<Bond>> get_bonds(std::unique_ptr<std::vector<Atom>> &atoms) {

    auto bonds = std::make_unique<std::vector<Bond>>();

    std::map<std::string, const Atom *> residue_atoms;

    auto current_residue_id = (*atoms)[0].residue_id();
    auto current_chain = (*atoms)[0].chain_id();

    for (const auto &atom: *atoms) {
        auto id = atom.residue_id();
        auto chain = atom.chain_id();

        /* Atom lies in the same residue */
        if (current_residue_id == id and current_chain == chain) {
            residue_atoms[atom.name()] = &atom;
        } else {
            /* New residue found, process the old one */
            update_bonds(bonds, residue_atoms);

            current_residue_id = id;
            current_chain = chain;
            residue_atoms.clear();
            residue_atoms[atom.name()] = &atom;
        }
    }

    update_bonds(bonds, residue_atoms);

    /* Add bonds on the protein backbone */
    std::map<int, const Atom *> N_backbone;
    std::map<int, const Atom *> C_backbone;

    for (const auto &atom: *atoms) {
        /* Skip HETATM */
        if (not atom.hetatm()) {
            if (atom.name() == "C") {
                C_backbone[atom.residue_id()] = &atom;
            } else if (atom.name() == "N") {
                N_backbone[atom.residue_id()] = &atom;
            }
        }
    }

    for (const auto &[residue_id, atom]: C_backbone) {
        auto it = N_backbone.find(residue_id + 1);
        if (it != N_backbone.end()) {
            bonds->emplace_back(atom, it->second, 1);
        }
    }

    return bonds;
}
