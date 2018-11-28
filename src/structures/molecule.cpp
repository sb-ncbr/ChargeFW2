//
// Created by krab1k on 23/10/18.
//

#include "atom.h"
#include "bond.h"
#include "molecule.h"
#include <utility>
#include <tuple>
#include <map>
#include <string>
#include <vector>
#include <iostream>


Molecule::Molecule(std::string name, std::unique_ptr<std::vector<Atom> > atoms,
                   std::unique_ptr<std::vector<Bond> > bonds, const std::map<int, int> &charges) {
    name_ = std::move(name);
    atoms_ = std::move(atoms);
    bonds_ = std::move(bonds);
    for (auto &[atom_no, charge]: charges) {
        (*atoms_)[atom_no].formal_charge_ = charge;
    }

    size_t n = atoms_->size();
    bond_info_.resize(n * n);
    for(const auto &bond: *bonds_) {
        int i = bond.first().index();
        int j = bond.second().index();
        char order = static_cast<char>(bond.order());
        bond_info_[i * n + j] = order;
        bond_info_[j * n + i] = order;
    }
}

std::ostream &operator<<(std::ostream &str, const Molecule &molecule) {
    return str << "Molecule: " << molecule.name_ << " Atoms: " << molecule.atoms_->size() << " Bonds: "
               << molecule.bonds_->size();
}

bool Molecule::bonded(const Atom &atom1, const Atom &atom2) const {
    size_t n = atoms_->size();
    return bond_info_[atom1.index() * n + atom2.index()] > 0;
}
