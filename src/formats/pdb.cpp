//
// Created by krab1k on 24.1.19.
//

#include <fstream>
#include <vector>
#include <fmt/format.h>
#include <boost/algorithm/string.hpp>

#include "config.h"
#include "pdb.h"
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

        /* Read header record and the pdbid */
        std::getline(file, line);
        std::string name = line.substr(62, 4);

        size_t idx = 0;
        while (std::getline(file, line)) {
            if (boost::starts_with(line, "ATOM")) {
                std::string atom_name = line.substr(12, 4);
                boost::trim(atom_name);
                std::string residue = line.substr(17, 3);
                size_t residue_id = std::stoul(line.substr(22, 4));
                double x = std::stod(line.substr(30, 8));
                double y = std::stod(line.substr(38, 8));
                double z = std::stod(line.substr(46, 8));
                std::string element_symbol = line.substr(76, 2);
                boost::trim(element_symbol);

                auto element = PeriodicTable::pte().getElement(element_symbol);

                atoms->emplace_back(idx, element, x, y, z, atom_name, residue_id, residue);
                idx++;
            }
        }

        auto bonds = std::make_unique<std::vector<Bond>>();
        std::map<size_t, int> charges;
        molecules->emplace_back(name, std::move(atoms), std::move(bonds), charges);

    } catch (const std::exception &) {
        fmt::print(stderr, "Invalid PDB/PQR file\n");
        exit(EXIT_FILE_ERROR);
    }

    return MoleculeSet(std::move(molecules));
}
