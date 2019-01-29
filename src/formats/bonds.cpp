//
// Created by krab1k on 28.1.19.
//

#include <vector>
#include <memory>
#include <string>
#include <tuple>
#include <fstream>
#include <map>
#include <sstream>

#include "bonds.h"
#include "../structures/atom.h"
#include "../structures/bond.h"
#include "config.h"


std::map<std::string, std::vector<std::tuple<std::string, std::string, int>>> load_residues_info();


void update_bonds(std::unique_ptr<std::vector<Bond>> &bonds, std::unique_ptr<std::vector<Atom>> &atoms,
                  std::map<std::string, const Atom *> residue_atoms);


std::map<std::string, std::vector<std::tuple<std::string, std::string, int>>> load_residues_info() {
    std::string filename(std::string(INSTALL_DIR) + "/share/amino_acids.txt");
    std::ifstream file(filename);
    if (!file) {
        fmt::print(stderr, "Unable to open amino acids data file: {}\n", filename);
        exit(EXIT_INTERNAL_ERROR);
    }

    std::map<std::string, std::vector<std::tuple<std::string, std::string, int>>> residues_info;

    std::string line;
    while (std::getline(file, line)) {
        std::string residue = line;
        while (not line.empty()) {
            std::string atom1_name;
            std::string atom2_name;
            int bond_order;

            std::getline(file, line);
            std::stringstream ss(line);
            ss >> atom1_name >> atom2_name >> bond_order;

            residues_info[residue].emplace_back(atom1_name, atom2_name, bond_order);
        }
    }
    return residues_info;
}


void update_bonds(std::unique_ptr<std::vector<Bond>> &bonds, std::unique_ptr<std::vector<Atom>> &atoms,
                  std::map<std::string, const Atom *> residue_atoms) {

    static auto residues_data = load_residues_info();

    auto residue = residue_atoms.begin()->second->residue();
    auto it = residues_data.find(residue);
    if (it != residues_data.end()) {
        for (const auto &[atom1_name, atom2_name, order]: it->second) {
            auto it1 = residue_atoms.find(atom1_name);
            auto it2 = residue_atoms.find(atom2_name);

            if (it1 != residue_atoms.end() and it2 != residue_atoms.end()) {
                auto atom1_idx = it1->second->index();
                auto atom2_idx = it2->second->index();
                bonds->emplace_back(&((*atoms)[atom1_idx]), &((*atoms)[atom2_idx]), order);
            }
        }
    }
}


std::unique_ptr<std::vector<Bond>> get_bonds(std::unique_ptr<std::vector<Atom>> &atoms) {

    auto bonds = std::make_unique<std::vector<Bond>>();

    std::map<std::string, const Atom *> residue_atoms;

    auto current_residue_id = (*atoms)[0].residue_id();

    for (const auto &atom: *atoms) {
        auto id = atom.residue_id();
        /* Atom lies in the same residue */
        if (current_residue_id == id) {
            residue_atoms[atom.name()] = &atom;
        } else {
            /* New residue found, process the old one */
            update_bonds(bonds, atoms, residue_atoms);

            current_residue_id = id;
            residue_atoms.clear();
            residue_atoms[atom.name()] = &atom;
        }
    }

    update_bonds(bonds, atoms, residue_atoms);

    /* Add bonds on the protein backbone */
    std::map<size_t, const Atom *> N_backbone;
    std::map<size_t, const Atom *> C_backbone;

    for (const auto &atom: *atoms) {
        if (atom.name() == "C") {
            C_backbone[atom.residue_id()] = &atom;
        } else if (atom.name() == "N") {
            N_backbone[atom.residue_id()] = &atom;
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
