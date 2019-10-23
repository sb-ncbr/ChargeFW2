//
// Created by krab1k on 24.1.19.
//

#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <set>
#include <fmt/format.h>
#include <boost/algorithm/string.hpp>

#include "common.h"
#include "chargefw2.h"
#include "../charges.h"
#include "../periodic_table.h"
#include "mol2.h"


MoleculeSet Mol2::read_file(const std::string &filename) {
    std::ifstream file(filename);
    if (!file) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    auto molecules = std::make_unique<std::vector<Molecule> >();

    std::string line;

    std::set<std::string> molecule_names;

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

            name = sanitize_name(name);
            name = get_unique_name(name, molecule_names);
            molecule_names.insert(name);

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

                size_t residue_id = 0;
                std::string residue = "UNL";
                as >> idx >> atom_name >> x >> y >> z >> atom_type;
                as >> residue_id >> residue;

                std::string element_symbol;
                auto it = atom_type.find('.');
                if (it != std::string::npos) {
                    element_symbol = atom_type.substr(0, it);
                } else {
                    element_symbol = atom_type;
                }

                auto element = PeriodicTable::pte().get_element_by_symbol(element_symbol);

                atoms->emplace_back(i, element, x, y, z, atom_name, residue_id, residue, "", false);
            }

            /* Read @<TRIPOS>BOND */
            std::getline(file, line);
            if (line != "@<TRIPOS>BOND") {
                fmt::print(stderr, "No BOND record\n");
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

            molecules->emplace_back(name, std::move(atoms), std::move(bonds));
        }
    }

    catch (const std::exception &) {
        fmt::print(stderr, "Invalid Mol2 file\n");
        exit(EXIT_FILE_ERROR);
    }

    return MoleculeSet(std::move(molecules));
}


void Mol2::save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    auto file = std::fopen(filename.c_str(), "w");
    if (!file) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    for (const auto &molecule: ms.molecules()) {
        try {
            auto chg = charges[molecule.name()];
            fmt::print(file, "@<TRIPOS>MOLECULE\n");
            fmt::print(file, "{}\n", molecule.name());
            fmt::print(file, "{} {}\n", molecule.atoms().size(), molecule.bonds().size());

            /* Try to guess if the molecule is protein or not */
            if (molecule.atoms()[0].chain_id().empty()) {
                fmt::print(file, "SMALL\n");
            } else {
                fmt::print(file, "PROTEIN\n");
            }
            fmt::print(file, "USER_CHARGES\n");
            fmt::print(file, "****\n");
            fmt::print(file, "Charges calculated by ChargeFW2 {}, method: {}\n", VERSION, charges.method_name());

            fmt::print(file, "@<TRIPOS>ATOM\n");
            for (size_t i = 0; i < molecule.atoms().size(); i++) {
                const auto &atom = molecule.atoms()[i];
                fmt::print(file, "{:>5d} {:<6s} {:>8.3f} {:>8.3f} {:>8.3f} {:s} {:>3d} {:>3s} {:>6.3f}\n",
                           i + 1, atom.name(), atom.pos()[0], atom.pos()[1], atom.pos()[2], atom.element().symbol(),
                           atom.residue_id(), atom.residue(), chg[i]);
            }

            fmt::print(file, "@<TRIPOS>BOND\n");
            for (size_t i = 0; i < molecule.bonds().size(); i++) {
                const auto &bond = molecule.bonds()[i];
                fmt::print(file, "{:>5d} {:>5d} {:>5d} {:>2d}\n",
                           i + 1, bond.first().index() + 1, bond.second().index() + 1, bond.order());
            }


        }
        catch (std::out_of_range &) {
            /* Do nothing */
        }
    }
    fclose(file);
}
