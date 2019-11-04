//
// Created by krab1k on 24/10/18.
//

#include <map>
#include <algorithm>
#include <memory>
#include <tuple>
#include <string>
#include <fmt/format.h>

#include "chargefw2.h"
#include "molecule_set.h"
#include "../periodic_table.h"
#include "../parameters.h"
#include "../method.h"


MoleculeSet::MoleculeSet(std::unique_ptr<std::vector<Molecule> > molecules) : molecules_{std::move(molecules)} {
    for (auto &molecule: *molecules_) {
        for (auto &atom: *molecule.atoms_)
            atom.molecule_ = &molecule;

        for (auto &bond: *molecule.bonds_)
            bond.molecule_ = &molecule;
    }
}


void MoleculeSet::info() const {
    fmt::print("Number of molecules: {}\n", molecules_->size());
    std::map<size_t, int> counts;
    size_t n_atoms = 0;
    for (const auto &m: *molecules_) {
        for (auto &a : m.atoms()) {
            counts[a.atom_type()] += 1;
            n_atoms++;
        }
    }

    fmt::print("Number of atoms: {}\n", n_atoms);
    for (auto &[key, val]: counts) {
        auto[symbol, cls, type] = atom_types_[key];
        fmt::print("{:2s} {} {}: {}\n", symbol, cls, type, val);
    }
}


void MoleculeSet::classify_atoms(AtomClassifier cls) {
    switch (cls) {
        case AtomClassifier::PLAIN: {
            for (auto &molecule: *molecules_) {
                for (auto &atom: *molecule.atoms_) {
                    auto tuple = std::make_tuple(atom.element().symbol(), std::string("plain"), std::string("*"));
                    set_atom_type(atom, tuple);
                }
            }
            break;
        }
        case AtomClassifier::HBO: {
            for (auto &molecule: *molecules_) {
                std::vector<int> max_bond_orders = get_max_bond_orders(molecule);
                for (auto &atom: *molecule.atoms_) {
                    auto tuple = std::make_tuple(atom.element().symbol(), "hbo",
                                                 std::to_string(max_bond_orders[atom.index_]));
                    set_atom_type(atom, tuple);
                }
            }
            break;
        }
    }
}


void MoleculeSet::classify_bonds(BondClassifier cls) {
    switch (cls) {
        case BondClassifier::PLAIN: {
            for (auto &molecule: *molecules_) {
                for (auto &bond: *molecule.bonds_) {
                    auto tuple = std::make_tuple(bond.first().element().symbol(), bond.second().element().symbol(),
                                                 std::string("plain"), std::string("*"));
                    set_bond_type(bond, tuple);
                }
            }
            break;
        }

        case BondClassifier::BO: {
            for (auto &molecule: *molecules_) {
                for (auto &bond: *molecule.bonds_) {
                    auto tuple = std::make_tuple(bond.first().element().symbol(), bond.second().element().symbol(),
                                                 std::string("bo"), std::to_string(bond.order_));
                    set_bond_type(bond, tuple);
                }
            }
            break;
        }
    }
}


void MoleculeSet::set_atom_type(Atom &atom, const std::tuple<std::string, std::string, std::string> &tuple) {
    auto it = find(atom_types_.begin(), atom_types_.end(), tuple);
    if (it == atom_types_.end()) {
        atom_types_.push_back(tuple);
        atom.atom_type_ = atom_types_.size() - 1;
    } else {
        atom.atom_type_ = static_cast<size_t>(distance(atom_types_.begin(), it));
    }
}


void MoleculeSet::set_bond_type(Bond &bond,
                                const std::tuple<std::string, std::string, std::string, std::string> &tuple) {
    auto it = find(bond_types_.begin(), bond_types_.end(), tuple);
    if (it == bond_types_.end()) {
        bond_types_.push_back(tuple);
        bond.bond_type_ = bond_types_.size() - 1;
    } else {
        bond.bond_type_ = static_cast<size_t>(distance(bond_types_.begin(), it));
    }
}


inline bool check_bond_symbols(const Bond &bond, const std::string &symbol1, const std::string &symbol2) {
    if (symbol1 != "*" and symbol2 != "*") {
        if ((bond.first().element().symbol() != symbol1 or bond.second().element().symbol() != symbol2) and
            (bond.second().element().symbol() != symbol1 or bond.first().element().symbol() != symbol2)) {
            return false;
        }
    }
    return true;
}


size_t MoleculeSet::classify_bonds_from_parameters(const Parameters &parameters, bool remove_unclassified,
                                                   bool permissive_types) {
    std::vector<int> unclassified;
    bond_types_ = parameters.bond()->keys();
    int m = 0;
    for (auto &molecule: *molecules_) {
        for (auto &bond: *molecule.bonds_) {
            bool found = false;
            for (size_t i = 0; i < bond_types_.size(); i++) {
                auto &[symbol1, symbol2, cls, type] = bond_types_[i];
                if (not check_bond_symbols(bond, symbol1, symbol2)) {
                    continue;
                }
                if (cls == "plain") {
                    bond.bond_type_ = i;
                    found = true;
                    break;
                } else if (cls == "bo") {
                    if (std::to_string(bond.order_) == type) {
                        bond.bond_type_ = i;
                        found = true;
                        break;
                    }
                } else {
                    fmt::print(stderr, "BondClassifier {} not found\n", cls);
                    exit(EXIT_INTERNAL_ERROR);
                }
            }
            if (!found) {
                bool is_unclassified = true;
                if (permissive_types) {
                    for (size_t i = 0; i < bond_types_.size(); i++) {
                        auto &[symbol1, symbol2, cls, type] = bond_types_[i];
                        if (not check_bond_symbols(bond, symbol1, symbol2)) {
                            continue;
                        }
                        if (cls == "bo") {
                            /* Try to match smaller bond order */
                            if (std::to_string(bond.order_ - 1) == type) {
                                bond.bond_type_ = i;
                                is_unclassified = false;
                                break;
                            }
                        }
                    }
                }
                if (is_unclassified) {
                    unclassified.push_back(m);
                    break;
                }
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


size_t MoleculeSet::classify_atoms_from_parameters(const Parameters &parameters, bool remove_unclassified,
                                                   bool permissive_types) {
    std::vector<int> unclassified;
    atom_types_ = parameters.atom()->keys();
    int m = 0;
    for (auto &molecule: *molecules_) {
        std::vector<int> max_bond_orders = get_max_bond_orders(molecule);
        for (auto &atom: *molecule.atoms_) {
            bool found = false;
            for (size_t i = 0; i < atom_types_.size(); i++) {
                auto &[symbol, cls, type] = atom_types_[i];
                if (atom.element().symbol() != symbol) {
                    continue;
                }
                if (cls == "plain") {
                    atom.atom_type_ = i;
                    found = true;
                    break;
                } else if (cls == "hbo") {
                    auto current_type = std::to_string(max_bond_orders[atom.index_]);
                    if (current_type == type) {
                        atom.atom_type_ = i;
                        found = true;
                        break;
                    }
                } else {
                    fmt::print(stderr, "AtomClassifier {} not found\n", cls);
                    exit(EXIT_INTERNAL_ERROR);
                }
            }
            if (!found) {
                bool is_unclassified = true;
                if (permissive_types) {
                    for (size_t i = 0; i < atom_types_.size(); i++) {
                        auto &[symbol, cls, type] = atom_types_[i];
                        if (atom.element().symbol() != symbol) {
                            continue;
                        }
                        if (cls == "hbo") {
                            std::string current_type;
                            /* Try to match similar bond order */
                            if (max_bond_orders[atom.index_] == 0) {
                                current_type = "1";
                            } else {
                                current_type = std::to_string(max_bond_orders[atom.index_] - 1);
                            }
                            if (current_type == type) {
                                atom.atom_type_ = i;
                                is_unclassified = false;
                                break;
                            }
                        }
                    }
                }
                if (is_unclassified) {
                    unclassified.push_back(m);
                    break;
                }
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


std::vector<int> MoleculeSet::get_max_bond_orders(const Molecule &molecule) const {
    const size_t n = molecule.atoms().size();
    std::vector<int> max_bond_orders(n, 0);

    for (const auto &bond: molecule.bonds()) {
        size_t i1 = bond.first_->index_;
        size_t i2 = bond.second_->index_;
        int bo = bond.order_;
        max_bond_orders[i1] = std::max(max_bond_orders[i1], bo);
        max_bond_orders[i2] = std::max(max_bond_orders[i2], bo);
    }
    return max_bond_orders;
}


size_t MoleculeSet::classify_set_from_parameters(const Parameters &parameters, bool remove_unclassified,
                                                 bool permissive_types) {
    size_t unclassified = 0;
    if (parameters.atom() != nullptr)
        unclassified += classify_atoms_from_parameters(parameters, remove_unclassified, permissive_types);

    if (parameters.bond() != nullptr)
        unclassified += classify_bonds_from_parameters(parameters, remove_unclassified, permissive_types);

    return unclassified;
}


void MoleculeSet::fulfill_requirements(const std::vector<RequiredFeatures> &features) {
    for (const auto req: features) {
        switch (req) {
            case RequiredFeatures::BOND_DISTANCES: {
                for (auto &molecule: *molecules_) {
                    molecule.init_bond_info();
                    molecule.init_bond_distances();
                }
                break;
            }

            case RequiredFeatures::BOND_INFO: {
                for (auto &molecule: *molecules_) {
                    molecule.init_bond_info();
                }
                break;
            }

            case RequiredFeatures::DISTANCE_TREE: {
                for (auto &molecule: *molecules_) {
                    molecule.init_distance_tree();
                }
                break;
            }
        }
    }
}


bool MoleculeSet::has_proteins() const {
    return std::any_of(molecules_->begin(), molecules_->end(), [](const Molecule &m) { return m.is_protein(); });
}
