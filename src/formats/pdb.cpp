//
// Created by krab1k on 24.1.19.
//

#include <fstream>
#include <vector>
#include <fmt/format.h>
#include <boost/algorithm/string.hpp>

#include "chargefw2.h"
#include "../config.h"
#include "pdb.h"
#include "common.h"
#include "bonds.h"
#include "../periodic_table.h"


MoleculeSet PDB::read_file(const std::string &filename) {
    std::ifstream file(filename);
    if (!file) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    auto molecules = std::make_unique<std::vector<Molecule> >();

    std::string line;
    try {

        auto atoms = std::make_unique<std::vector<Atom>>();

        std::string name = filename;

        size_t idx = 0;
        bool atom_block_found = false;
        while (std::getline(file, line)) {
            if (boost::starts_with(line, "ATOM") or (config::read_hetatm and boost::starts_with(line, "HETATM"))) {
                atom_block_found = true;

                std::string atom_name = line.substr(12, 4);
                boost::trim(atom_name);
                auto residue = line.substr(17, 3);
                if (config::ignore_water and residue == "HOH") {
                    continue;
                }

                bool hetatm = false;
                if (line[0] == 'H') {
                    hetatm = true;
                }

                auto residue_id = std::stoi(line.substr(22, 4));
                auto chain_id = line.substr(21, 1);
                auto x = std::stod(line.substr(30, 8));
                auto y = std::stod(line.substr(38, 8));
                auto z = std::stod(line.substr(46, 8));
                auto symbol = line.substr(76, 2);

                auto element = PeriodicTable::pte().getElement(get_element_symbol(symbol));

                atoms->emplace_back(idx, element, x, y, z, atom_name, residue_id, residue, chain_id, hetatm);
                idx++;
            } else if (boost::starts_with(line, "TER")) {
                /* Chain ended, look for next one */
                continue;
            } else if (atom_block_found) {
                /* We are done reading ATOM and HETATM records */
                break;
            }
        }

        auto bonds = get_bonds(atoms);
        std::map<size_t, int> charges;
        molecules->emplace_back(name, std::move(atoms), std::move(bonds), charges);

    } catch (const std::exception &) {
        fmt::print(stderr, "Invalid PDB file\n");
        exit(EXIT_FILE_ERROR);
    }

    return MoleculeSet(std::move(molecules));
}
