//
// Created by krab1k on 24/10/18.
//

#include <map>
#include <algorithm>
#include <iostream>
#include <memory>
#include <tuple>
#include <string>

#include "molecule_set.h"
#include "../classifier.h"
#include "../periodic_table.h"
#include "../parameters.h"
#include "config.h"


MoleculeSet::MoleculeSet(std::unique_ptr<std::vector<Molecule> > molecules) : molecules_{std::move(molecules)} {
    for (auto &molecule: *molecules_) {
        for (auto &atom: *molecule.atoms_)
            atom.molecule_ = &molecule;

        for (auto &bond: *molecule.bonds_)
            bond.molecule_ = &molecule;
    }
}


void MoleculeSet::info() const {
    std::map<size_t, int> counts;
    for (const Molecule &m: *molecules_) {
        for (auto &a : m.atoms()) {
            counts[a.atom_type()] += 1;
        }
    }

    for (auto &[key, val]: counts) {
        auto[symbol, cls, type] = atom_types_[key];
        std::cout << symbol << " " << cls << " " << type << ": "
                  << val << std::endl;

    }
}


void MoleculeSet::classify_atoms(const Classifier &cls) {
    for (auto &molecule: *molecules_) {
        for (auto &atom: *molecule.atoms_) {
            auto type = cls.get_type(atom);
            auto tuple = std::make_tuple(atom.element().symbol(), cls.name(), type);
            auto it = std::find(atom_types_.begin(), atom_types_.end(), tuple);
            if (it == atom_types_.end()) {
                atom_types_.push_back(tuple);
                atom.atom_type_ = atom_types_.size() - 1;
            } else {
                atom.atom_type_ = static_cast<size_t>(std::distance(atom_types_.begin(), it));
            }
        }
    }
}


void MoleculeSet::classify_atoms_from_parameters(const Parameters &parameters) {
    atom_types_ = parameters.atom()->keys();
    for (auto &molecule: *molecules_) {
        for (auto &atom: *molecule.atoms_) {
            bool found = false;
            for (size_t i = 0; i < atom_types_.size(); i++) {
                auto &[symbol, cls, type] = atom_types_[i];
                if (atom.element().symbol() != symbol)
                    continue;
                if (cls == "plain") {
                    atom.atom_type_ = i;
                    found = true;
                    break;
                } else if (cls == "hbo") {
                    auto hbo = HBOClassifier();
                    auto current_type = hbo.get_type(atom);
                    if (current_type == type) {
                        atom.atom_type_ = i;
                        found = true;
                        break;
                    }
                } else {
                    std::cerr << "Classifier " << cls << " not found" << std::endl;
                    exit(EXIT_INTERNAL_ERROR);
                }
            }
            if (!found) {
                std::cerr << "No parameters for atom " << atom.element().symbol() << " in molecule "
                          << molecule.name() << std::endl;
                exit(EXIT_GENERAL_ERROR);
            }
        }
    }
}
