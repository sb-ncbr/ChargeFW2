//
// Created by krab1k on 24/10/18.
//

#include <string>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>
#include <memory>
#include <boost/algorithm/string.hpp>

#include "../structures/Atom.h"
#include "../structures/Bond.h"
#include "../structures/Molecule.h"
#include "../PeriodicTable.h"
#include "SDF.h"


MoleculeSet SDF::read_file(const std::string &filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
	
    auto molecules = std::make_unique<std::vector<Molecule> >();
    try {
        while (std::getline(file, line)) {
            std::string name = line; // Read name of the molecule
            std::getline(file, line); // Line with comments
            std::getline(file, line); // Line with comments

            std::getline(file, line); // Line with counts

            size_t n_atoms = std::stoul(line.substr(0, 3));
            size_t n_bonds = std::stoul(line.substr(3, 3));

            auto atoms = std::make_unique<std::vector<Atom> >();
            atoms->reserve(n_atoms);

            for (unsigned i = 0; i < n_atoms; i++) {
                std::getline(file, line);
                double x = std::stod(line.substr(0, 10));
                double y = std::stod(line.substr(10, 10));
                double z = std::stod(line.substr(20, 10));

                auto element_symbol = boost::trim_copy(line.substr(31, 3));
                boost::to_lower(element_symbol);
                element_symbol[0] = static_cast<char>(std::toupper(element_symbol[0]));

                auto element = PeriodicTable::pte().getElement(element_symbol);

                atoms->emplace_back(i, element, x, y, z);
            }

            auto bonds = std::make_unique<std::vector<Bond> >();
            bonds->reserve(n_bonds);

            for (unsigned i = 0; i < n_bonds; i++) {
                std::getline(file, line);
                int first = std::stoi(line.substr(0, 3));
                int second = std::stoi(line.substr(3, 3));
                int order = std::stoi(line.substr(6, 3));

                bonds->emplace_back(&((*atoms)[first - 1]), &((*atoms)[second - 1]), order);
            }

            std::map<int, int> charges;
            do {
                std::getline(file, line);
                if(line.substr(0, 6) == "M  CHG") {
                    int count = std::stoi(line.substr(6, 3));
                    const int base = 9;
                    for(int i = 0; i < count; i++) {
                        int atom_no = std::stoi(line.substr(base + i * 8, 4));
                        int charge = std::stoi(line.substr(base + i * 8 + 4, 4));
                        charges[atom_no - 1] = charge;
                    }
                }
            } while (line != "$$$$");

            molecules->emplace_back(name, std::move(atoms), std::move(bonds), charges);
        }
    }
    catch (const std::invalid_argument &e) {
        std::cerr << "Invalid SDF file" << std::endl;
        exit(EXIT_FAILURE);
    }
    return MoleculeSet(std::move(molecules));
}
