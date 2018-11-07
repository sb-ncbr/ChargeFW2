//
// Created by krab1k on 24/10/18.
//

#include <map>
#include <iostream>
#include <memory>
#include <tuple>
#include <string>

#include "MoleculeSet.h"
#include "../Classifier.h"
#include "../PeriodicTable.h"
#include "../Parameters.h"

MoleculeSet::MoleculeSet(std::unique_ptr<std::vector<Molecule> > molecules) : molecules_{std::move(molecules)} {
    for (auto &molecule: *molecules_) {
        for (auto &atom: *molecule.atoms_)
            atom.molecule_ = &molecule;

        for (auto &bond: *molecule.bonds_)
            bond.molecule_ = &molecule;
    }
}


void MoleculeSet::info() const {
    std::map<std::tuple<std::string, std::string, std::string>, int> counts;
    for (const Molecule &m: *molecules_) {
        for (auto &a : m.atoms()) {
            counts[a.atom_type()] += 1;
        }
    }

    for (auto &[key, val]: counts) {
        auto[symbol, cls, type] = key;
        std::cout << symbol << " " << cls << " " << type << ": "
                  << val << std::endl;

    }
}


void MoleculeSet::classify_atoms(std::string classifier) {
    auto cls = std::unique_ptr<Classifier>(nullptr);

    if (classifier == "plain") {
        cls = std::make_unique<PlainClassifier>();
    } else if (classifier == "hbo") {
        cls = std::make_unique<HBOClassifier>();
    } else {
        std::cerr << "Unknown classifier " << classifier << std::endl;
        exit(EXIT_FAILURE);
    }

    for (auto &molecule: *molecules_) {
        for (auto &atom: *molecule.atoms_) {
            atom.atom_type_ = std::make_tuple(atom.element().symbol(), cls->name(), cls->get_type(atom));
        }
    }
}


void MoleculeSet::classify_atoms_from_parameters(const Parameters &parameters) {
    for (auto &molecule: *molecules_) {
        for (auto &atom: *molecule.atoms_) {
            bool found = false;
            for (const auto &key: parameters.atom()->keys()) {
                auto &[symbol, cls, type] = key;
                if (atom.element().symbol() != symbol)
                    continue;
                if (cls == "plain") {
                    atom.atom_type_ = std::make_tuple(atom.element().symbol(), "plain", "*");
                    found = true;
                    break;
                } else if (cls == "hbo") {
                    auto hbo = HBOClassifier();
                    auto current_type = hbo.get_type(atom);
                    if (current_type == type) {
                        atom.atom_type_ = std::make_tuple(atom.element().symbol(), "hbo", current_type);
                        found = true;
                        break;
                    }
                } else {
                    std::cerr << "Classifier " << cls << " not found" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            if (!found) {
                std::cerr << "No parameters for atom " << atom.element().symbol() << " in molecule "
                          << molecule.name() << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
}
