//
// Created by krab1k on 24.1.19.
//

#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <fmt/format.h>
#include <boost/algorithm/string.hpp>

#include "../periodic_table.h"
#include "mol2.h"
#include "config.h"


MoleculeSet Mol2::read_file(const std::string &filename) {
    std::ifstream file(filename);
    if (!file) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    auto molecules = std::make_unique<std::vector<Molecule> >();

    std::string line;

    try {
        while (std::getline(file, line)) {

            /* Skip comments or empty lines*/
            while (boost::starts_with(line, "#") || line.empty() || line != "@<TRIPOS>MOLECULE") {
                std::getline(file, line);
            }

            if (line != "@<TRIPOS>MOLECULE") {
                fmt::print(stderr, "No MOLECULE record\n");
                throw std::exception();
            }

            std::string name;
            std::getline(file, name);

            std::getline(file, line);

            std::stringstream ss(line);
            size_t n_atoms, n_bonds;
            ss >> n_atoms >> n_bonds;

            /* Skip the rest until atom records */
            do {
                std::getline(file, line);
            } while (line != "@<TRIPOS>ATOM");

            auto atoms = std::make_unique<std::vector<Atom> >();
            atoms->reserve(n_atoms);

            for (size_t i = 0; i < n_atoms; i++) {
                std::getline(file, line);
                std::stringstream as(line);
                size_t idx;
                std::string atom_name, atom_type;
                double x, y, z;
                as >> idx >> atom_name >> x >> y >> z >> atom_type;

                std::string element_symbol;
                auto it = atom_type.find('.');
                if (it != std::string::npos) {
                    element_symbol = atom_type.substr(0, it);
                } else {
                    element_symbol = atom_type;
                }

                auto element = PeriodicTable::pte().getElement(element_symbol);

                atoms->emplace_back(i, element, x, y, z);
            }

            /* Read @<TRIPOS>BOND */
            std::getline(file, line);
            if (line != "@<TRIPOS>BOND") {
                fmt::print(stderr, "No MOLECULE record\n");
                throw std::exception();
            }

            auto bonds = std::make_unique<std::vector<Bond> >();
            bonds->reserve(n_bonds);

            for (size_t i = 0; i < n_bonds; i++) {
                std::getline(file, line);
                std::stringstream bs(line);

                size_t idx, first, second;
                std::string type;
                bs >> idx >> first >> second >> type;

                int order;
                try {
                    order = std::stoi(type);
                } catch (std::invalid_argument &) {
                    /* Set bond order to 1 for aromatic and other bond types */
                    order = 1;
                }
                bonds->emplace_back(&((*atoms)[first - 1]), &((*atoms)[second - 1]), order);
            }

            std::map<size_t, int> charges;
            molecules->emplace_back(name, std::move(atoms), std::move(bonds), charges);
        }
    }

    catch (const std::exception &e) {
        fmt::print(stderr, "Invalid Mol2 file\n");
        exit(EXIT_FILE_ERROR);
    }

    return MoleculeSet(std::move(molecules));
}
