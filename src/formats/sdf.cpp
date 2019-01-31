//
// Created by krab1k on 24/10/18.
//

#include <string>
#include <fstream>
#include <sstream>
#include <tuple>
#include <vector>
#include <memory>
#include <fmt/format.h>
#include <boost/algorithm/string.hpp>

#include "../structures/atom.h"
#include "../structures/bond.h"
#include "../structures/molecule.h"
#include "../periodic_table.h"
#include "sdf.h"
#include "config.h"


MoleculeSet SDF::read_file(const std::string &filename) {
    std::ifstream file(filename);
    if (!file) {
        fmt::print("Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    std::string line;
	
    auto molecules = std::make_unique<std::vector<Molecule> >();
    try {
        while (std::getline(file, line)) {
            std::string name = line; // Read name of the molecule
            std::getline(file, line); // Line with comments
            std::getline(file, line); // Line with comments

            std::getline(file, line); // Line with counts

            std::string version = line.substr(34, 5);
            if (version == "V2000") {
                size_t n_atoms = std::stoul(line.substr(0, 3));
                size_t n_bonds = std::stoul(line.substr(3, 3));

                auto atoms = std::make_unique<std::vector<Atom> >();
                atoms->reserve(n_atoms);

                for (size_t i = 0; i < n_atoms; i++) {
                    std::getline(file, line);
                    double x = std::stod(line.substr(0, 10));
                    double y = std::stod(line.substr(10, 10));
                    double z = std::stod(line.substr(20, 10));

                    auto element_symbol = boost::trim_copy(line.substr(31, 3));
                    boost::to_lower(element_symbol);
                    element_symbol[0] = static_cast<char>(std::toupper(element_symbol[0]));

                    auto element = PeriodicTable::pte().getElement(element_symbol);

                    atoms->emplace_back(i, element, x, y, z, element->symbol(), 0, "UNL", "");
                }

                auto bonds = std::make_unique<std::vector<Bond> >();
                bonds->reserve(n_bonds);

                for (size_t i = 0; i < n_bonds; i++) {
                    std::getline(file, line);
                    size_t first = std::stoul(line.substr(0, 3));
                    size_t second = std::stoul(line.substr(3, 3));
                    int order = std::stoi(line.substr(6, 3));

                    bonds->emplace_back(&((*atoms)[first - 1]), &((*atoms)[second - 1]), order);
                }

                std::map<size_t, int> charges;
                do {
                    std::getline(file, line);
                    if(line.substr(0, 6) == "M  CHG") {
                        size_t count = std::stoul(line.substr(6, 3));
                        const size_t base = 9;
                        for(size_t i = 0; i < count; i++) {
                            size_t atom_no = std::stoul(line.substr(base + i * 8, 4));
                            int charge = std::stoi(line.substr(base + i * 8 + 4, 4));
                            charges[atom_no - 1] = charge;
                        }
                    }
                } while (line != "$$$$");

                molecules->emplace_back(name, std::move(atoms), std::move(bonds), charges);
            }
            else if (version == "V3000") {
                /* Skip 'M  V30 BEGIN CTAB' line */
                std::getline(file, line);

                /* Read atom & bond counts */
                std::getline(file, line);
                size_t n_atoms;
                size_t n_bonds;
                std::stringstream ss(line.substr(14));
                ss >> n_atoms >> n_bonds;

                auto atoms = std::make_unique<std::vector<Atom> >();
                atoms->reserve(n_atoms);

                std::map<size_t, int> charges;

                /* Skip 'M  V30 BEGIN ATOM' line */
                std::getline(file, line);

                /* Read info about the atoms */
                for (size_t i = 0; i < n_atoms; i++) {
                    std::getline(file, line);

                    ss.str(line.substr(7));
                    ss.clear();

                    std::string element_symbol;
                    size_t idx;
                    double x, y, z;
                    int formal_charge = 0;
                    ss >> idx >> element_symbol >> x >> y >> z;

                    /* Search for charge data */
                    auto pos = line.find("CHG=");
                    if (pos != std::string::npos) {
                        ss.str(line.substr(pos + 4));
                        ss.clear();
                        ss >> formal_charge;
                    }

                    /* Should we continue to next line */
                    while (line.back() == '-') {
                        std::getline(file, line);
                        pos = line.find("CHG=");
                        if (pos != std::string::npos) {
                            ss.str(line.substr(pos + 4));
                            ss.clear();
                            ss >> formal_charge;
                        }
                    }

                    charges[idx - 1] = formal_charge;

                    boost::to_lower(element_symbol);
                    element_symbol[0] = static_cast<char>(std::toupper(element_symbol[0]));

                    auto element = PeriodicTable::pte().getElement(element_symbol);

                    atoms->emplace_back(i, element, x, y, z, element->symbol(), 0, "UNL", "");
                }

                /* Skip 'M  V30 END ATOM' line */
                std::getline(file, line);

                /* Skip 'M  V30 BEGIN BOND' line */
                std::getline(file, line);

                auto bonds = std::make_unique<std::vector<Bond> >();
                bonds->reserve(n_bonds);

                for(size_t i = 0; i < n_bonds; i++) {
                    std::getline(file, line);
                    ss.str(line.substr(7));
                    ss.clear();
                    size_t idx, first, second;
                    int order;
                    ss >> idx >> order >> first >> second;

                    bonds->emplace_back(&((*atoms)[first - 1]), &((*atoms)[second - 1]), order);
                }

                /* Skip everything until the end of the record */
                do {
                    std::getline(file, line);
                } while (line != "$$$$");

                molecules->emplace_back(name, std::move(atoms), std::move(bonds), charges);

            }
            else {
                fmt::print(stderr, "Invalid MOL version {} inside SDF file\n", version);
                exit(EXIT_FILE_ERROR);
            }
        }
    }
    catch (const std::exception &) {
        fmt::print(stderr, "Invalid SDF file\n");
        exit(EXIT_FILE_ERROR);
    }
    return MoleculeSet(std::move(molecules));
}
