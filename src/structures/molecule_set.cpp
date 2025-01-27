#include <map>
#include <algorithm>
#include <memory>
#include <tuple>
#include <string>
#include <type_traits>
#include <fmt/format.h>

#include "chargefw2.h"
#include "molecule.h"
#include "molecule_set.h"
#include "../parameters.h"
#include "../method.h"
#include "../exceptions/internal_exception.h"

bool check_atom_type(const Molecule &molecule, const Atom &atom, const std::string &symbol, const std::string &cls,
                     const std::string &type, bool permissive);

bool check_bond_type(const Bond &bond, const std::string &cls, const std::string &type, bool permissive);

bool check_bond_atoms(const Molecule &molecule, const Bond &bond, const bond_t &bond_type);


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
            counts[a.type()] += 1;
            n_atoms++;
        }
    }

    fmt::print("Number of atoms: {}\n", n_atoms);
    for (auto &[key, val]: counts) {
        auto[symbol, cls, type] = atom_types_[key];
        fmt::print("{:2s} {:6s} {:4s}: {}\n", symbol, cls, type, val);
    }
}


void MoleculeSet::classify_atoms(AtomClassifier cls) {
    switch (cls) {
        case AtomClassifier::PLAIN: {
            for (auto &molecule: *molecules_) {
                for (auto &atom: *molecule.atoms_) {
                    auto tuple = std::make_tuple(atom.element().symbol(), std::string("plain"), std::string("*"));
                    set_type<Atom>(atom, tuple);
                }
            }
            break;
        }
        case AtomClassifier::HBO: {
            for (auto &molecule: *molecules_) {
                for (auto &atom: *molecule.atoms_) {
                    auto tuple = std::make_tuple(atom.element().symbol(), std::string("hbo"),
                                                 std::to_string(molecule.get_max_bond_orders()[atom.index_]));
                    set_type<Atom>(atom, tuple);
                }
            }
            break;
        }
        case AtomClassifier::BONDED: {
            for (auto &molecule: *molecules_) {
                for (size_t i = 0; i < molecule.atoms().size(); i++) {
                    auto &atom = (*molecule.atoms_)[i];
                    auto atom_type = molecule.get_bonded_elements()[i];
                    auto tuple = std::make_tuple(atom.element().symbol(), std::string("bonded"), atom_type);
                    set_type<Atom>(atom, tuple);
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
                    auto tuple = std::make_tuple(bond.first().element().symbol(), std::string("plain"), std::string("*"),
                                                 bond.second().element().symbol(), std::string("plain"), std::string("*"),
                                                 std::string("plain"), std::string("*"));
                    set_type<Bond>(bond, tuple);
                }
            }
            break;
        }

        case BondClassifier::BO: {
            for (auto &molecule: *molecules_) {
                for (auto &bond: *molecule.bonds_) {
                    auto tuple = std::make_tuple(bond.first().element().symbol(), std::string("plain"), std::string("*"),
                                                 bond.second().element().symbol(), std::string("plain"), std::string("*"),
                                                 std::string("bo"), std::to_string(bond.order_));
                    set_type<Bond>(bond, tuple);
                }
            }
            break;
        }
    }
}


bool check_atom_type(const Molecule &molecule, const Atom &atom, const std::string &symbol, const std::string &cls,
                     const std::string &type, bool permissive = false) {
    if (symbol != "*" and symbol != atom.element().symbol()) {
        return false;
    }

    std::string current_type;
    if (cls == "plain") {
        current_type = "*";
    } else if (cls == "hbo") {
        auto actual_order = molecule.get_max_bond_orders()[atom.index()];
        if (permissive) {
            current_type = actual_order == 0 ? std::string("1") : std::to_string(actual_order - 1);
        } else {
            current_type = std::to_string(actual_order);
        }
    } else if (cls == "bonded") {
        current_type = molecule.get_bonded_elements()[atom.index()];
    } else {
        throw InternalException(fmt::format("AtomClassifier {} not found", cls));
    }

    return current_type == type;
}


bool check_bond_type(const Bond &bond, const std::string &cls, const std::string &type, bool permissive = false) {
    std::string current_type;
    if (cls == "plain") {
        current_type = "*";
    } else if (cls == "bo") {
        if (permissive) {
            /* Try to match smaller bond order */
            current_type = std::to_string(bond.order() - 1);
        } else {
            current_type = std::to_string(bond.order());
        }

    } else {
        throw InternalException(fmt::format("BondClassifier {} not found", cls));
    }
    return current_type == type;
}


bool check_bond_atoms(const Molecule &molecule, const Bond &bond, const bond_t &bond_type) {

    auto check_bond_atoms_in_order = [&molecule, &bond_type](const Atom &atom1, const Atom &atom2) {
        const auto &[symbol1, cls1, type1, symbol2, cls2, type2, cls_b, type_b] = bond_type;
        return check_atom_type(molecule, atom1, symbol1, cls1, type1) and
               check_atom_type(molecule, atom2, symbol2, cls2, type2);
    };

    const auto &atom1 = bond.first();
    const auto &atom2 = bond.second();

    return check_bond_atoms_in_order(atom1, atom2) or check_bond_atoms_in_order(atom2, atom1);
}


template<typename AB, typename AB_t>
void MoleculeSet::set_type(AB &object, const AB_t &type) {

    std::vector<AB_t> *types;
    if constexpr (std::is_same_v<AB, Atom>) {
        types = &atom_types_;
    } else {
        types = &bond_types_;
    }

    auto it = find(types->begin(), types->end(), type);
    if (it == types->end()) {
        types->push_back(type);
        object.type_ = types->size() - 1;
    } else {
        object.type_ = static_cast<size_t>(distance(types->begin(), it));
    }
}


template<typename AB, typename AB_t>
size_t MoleculeSet::classify_objects_from_parameters(const Parameters &parameters,
                                                     bool remove_unclassified,
                                                     bool permissive_types) {

    std::vector<AB_t> *types;
    if constexpr (std::is_same_v<AB, Atom>) {
        atom_types_ = parameters.atom()->keys();
        types = &atom_types_;
    } else {
        bond_types_ = parameters.bond()->keys();
        types = &bond_types_;
    }

    std::vector<int> unclassified;
    int m = 0;
    for (auto &molecule: *molecules_) {

        std::vector<AB> *objects;
        if constexpr (std::is_same_v<AB, Atom>) {
            objects = molecule.atoms_.get();
        } else {
            objects = molecule.bonds_.get();
        }

        for (auto &object: *objects) {
            bool found = false;
            for (size_t i = 0; i < types->size(); i++) {
                bool cond;
                if constexpr (std::is_same_v<AB, Atom>) {
                    const auto &[symbol, cls, type] = (*types)[i];
                    cond = check_atom_type(molecule, object, symbol, cls, type);
                } else {
                    const auto bond_type = (*types)[i];
                    const auto &[symbol1, cls1, type1, symbol2, cls2, type2, cls_b, type_b] = bond_type;
                    if (not check_bond_atoms(molecule, object, bond_type)) {
                        continue;
                    }
                    cond = check_bond_type(object, cls_b, type_b);
                }
                if (cond) {
                    object.type_ = i;
                    found = true;
                    break;
                }
            }

            if (!found) {
                bool is_unclassified = true;
                if (permissive_types) {
                    for (size_t i = 0; i < types->size(); i++) {
                        bool cond;
                        if constexpr (std::is_same_v<AB, Atom>) {
                            const auto &[symbol, cls, type] = (*types)[i];
                            cond = check_atom_type(molecule, object, symbol, cls, type, true);
                        } else {
                            const auto bond_type = (*types)[i];
                            const auto &[symbol1, cls1, type1, symbol2, cls2, type2, cls_b, type_b] = bond_type;
                            if (not check_bond_atoms(molecule, object, bond_type)) {
                                continue;
                            }
                            cond = check_bond_type(object, cls_b, type_b, true);
                        }

                        if (cond) {
                            object.type_ = i;
                            is_unclassified = false;
                            break;
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


size_t MoleculeSet::classify_set_from_parameters(const Parameters &parameters, bool remove_unclassified,
                                                 bool permissive_types) {
    size_t unclassified = 0;
    if (parameters.atom() != nullptr) {
        unclassified += classify_objects_from_parameters<Atom>(parameters, remove_unclassified, permissive_types);
    }

    if (parameters.bond() != nullptr) {
        unclassified += classify_objects_from_parameters<Bond>(parameters, remove_unclassified, permissive_types);
    }

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
