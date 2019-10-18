//
// Created by krab1k on 24.1.19.
//

#include <fstream>
#include <vector>
#include <fmt/format.h>
#include <filesystem>
#include <boost/algorithm/string.hpp>

#include "chargefw2.h"
#include "../config.h"
#include "pdb.h"
#include "common.h"
#include "bonds.h"
#include "../periodic_table.h"

namespace fs = std::filesystem;


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

        std::string name = fs::path(filename).filename().replace_extension();

        size_t idx = 0;
        while (std::getline(file, line)) {
            if (boost::starts_with(line, "HEADER")) {
                name = line.substr(62, 4);
                continue;
            }

            if (boost::starts_with(line, "ATOM") or (config::read_hetatm and boost::starts_with(line, "HETATM"))) {
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

                auto alt_loc = line[16];
                auto element = PeriodicTable::pte().get_element_by_symbol(get_element_symbol(symbol));

                if (alt_loc == ' ' or not is_already_loaded(*atoms, atom_name, residue_id)) {
                    atoms->emplace_back(idx, element, x, y, z, atom_name, residue_id, residue, chain_id, hetatm);
                    idx++;

                }
            }

            /* If they are multiple models present, end after the first one */
            if (boost::starts_with(line, "ENDMDL")) {
                break;
            }
        }

        auto bonds = get_bonds(atoms);
        std::map<size_t, int> charges;
        molecules->emplace_back(sanitize_name(name), std::move(atoms), std::move(bonds), charges);

    } catch (const std::exception &) {
        fmt::print(stderr, "Invalid PDB file\n");
        exit(EXIT_FILE_ERROR);
    }

    return MoleculeSet(std::move(molecules));
}
