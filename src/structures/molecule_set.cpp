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
    std::cout << "Number of molecules: " << molecules_->size() << std::endl;
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


void MoleculeSet::classify_atoms(const AtomClassifier &cls) {
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


void MoleculeSet::classify_bonds(const BondClassifier &cls) {
    for (auto &molecule: *molecules_) {
        for (auto &bond: *molecule.bonds_) {
            auto type = cls.get_type(bond);
            auto tuple = std::make_tuple(bond.first().element().symbol(), bond.second().element().symbol(), cls.name(),
                                         type);
            auto it = std::find(bond_types_.begin(), bond_types_.end(), tuple);
            if (it == bond_types_.end()) {
                bond_types_.push_back(tuple);
                bond.bond_type_ = bond_types_.size() - 1;
            } else {
                bond.bond_type_ = static_cast<size_t>(std::distance(bond_types_.begin(), it));
            }
        }
    }
}


size_t MoleculeSet::classify_bonds_from_parameters(const Parameters &parameters, bool remove_unclassified) {
    std::vector<int> unclassified;
    bond_types_ = parameters.bond()->keys();
    int m = 0;
    for (auto &molecule: *molecules_) {
        for (auto &bond: *molecule.bonds_) {
            bool found = false;
            for (size_t i = 0; i < bond_types_.size(); i++) {
                auto &[symbol1, symbol2, cls, type] = bond_types_[i];
                if (symbol1 != "*" or symbol2 != "*") {
                    if ((bond.first().element().symbol() != symbol1 or bond.second().element().symbol() != symbol2) and
                        (bond.second().element().symbol() != symbol1 or bond.first().element().symbol() != symbol2)) {
                        continue;
                    }
                }
                if (cls == "plain") {
                    bond.bond_type_ = i;
                    found = true;
                    break;
                } else if (cls == "bo") {
                    auto bo = BOBondClassifier();
                    auto current_type = bo.get_type(bond);
                    if (current_type == type) {
                        bond.bond_type_= i;
                        found = true;
                        break;
                    }
                }
                else {
                    std::cerr << "BondClassifier " << cls << " not found" << std::endl;
                    exit(EXIT_INTERNAL_ERROR);
                }
            }
            if (!found) {
                unclassified.push_back(m);
                break;
            }
        }
        m++;
    }

    if (remove_unclassified) {
        // Need to iterate in reverse order to maintain indices correctness
        for (size_t i = 0; i < unclassified.size(); i++) {
            molecules_->erase(molecules_->begin() + unclassified[unclassified.size() - i - 1]);
        }
    }
    return unclassified.size();
}


size_t MoleculeSet::classify_atoms_from_parameters(const Parameters &parameters, bool remove_unclassified) {
    std::vector<int> unclassified;
    atom_types_ = parameters.atom()->keys();
    int m = 0;
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
                    auto hbo = HBOAtomClassifier();
                    auto current_type = hbo.get_type(atom);
                    if (current_type == type) {
                        atom.atom_type_ = i;
                        found = true;
                        break;
                    }
                } else {
                    std::cerr << "AtomClassifier " << cls << " not found" << std::endl;
                    exit(EXIT_INTERNAL_ERROR);
                }
            }
            if (!found) {
                unclassified.push_back(m);
                break;
            }
        }
        m++;
    }

    if (remove_unclassified) {
        // Need to iterate in reverse order to maintain indices correctness
        for (size_t i = 0; i < unclassified.size(); i++) {
            molecules_->erase(molecules_->begin() + unclassified[unclassified.size() - i - 1]);
        }
    }
    return unclassified.size();
}


size_t MoleculeSet::classify_set_from_parameters(const Parameters &parameters, bool remove_unclassified) {
    size_t unclassified = 0;
    if (parameters.atom() != nullptr)
        unclassified += classify_atoms_from_parameters(parameters, remove_unclassified);

    if (parameters.bond() != nullptr)
        unclassified += classify_bonds_from_parameters(parameters, remove_unclassified);

    return unclassified;
}
