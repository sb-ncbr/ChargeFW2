//
// Created by krab1k on 23/10/18.
//

#include "atom.h"
#include "bond.h"
#include "molecule.h"
#include "../geometry.h"

#include <queue>
#include <utility>
#include <tuple>
#include <map>
#include <string>
#include <vector>


Molecule::Molecule(std::string name, std::unique_ptr<std::vector<Atom> > atoms,
                   std::unique_ptr<std::vector<Bond> > bonds, const std::map<int, int> &charges) {
    name_ = std::move(name);
    atoms_ = std::move(atoms);
    bonds_ = std::move(bonds);
    for (auto &[atom_no, charge]: charges) {
        (*atoms_)[atom_no].formal_charge_ = charge;
    }
}


bool Molecule::bonded(const Atom &atom1, const Atom &atom2) const {
    const size_t n = atoms_->size();
    return bond_info_[atom1.index() * n + atom2.index()] > 0;
}


int Molecule::bond_order(const Atom &atom1, const Atom &atom2) const {
    const size_t n = atoms_->size();
    return static_cast<int>(bond_info_[atom1.index() * n + atom2.index()]);
}


int Molecule::degree(const Atom &atom) const {
    const size_t n = atoms_->size();
    int sum = 0;
    for (size_t i = 0; i < n; i++) {
        sum += bond_info_[atom.index() * n + i];
    }
    return sum;
}


std::vector<int> Molecule::get_bonded(int atom_idx) const {
    const size_t n = atoms_->size();
    std::vector<int> res;

    for (size_t j = 0; j < n; j++) {
        if (bond_info_[atom_idx * n + j]) {
            res.push_back(static_cast<int>(j));
        }
    }
    return res;
}


void Molecule::init_bond_info() {
    const size_t n = atoms_->size();
    bond_info_.resize(n * n);

    for (const auto &bond: *bonds_) {
        int i = bond.first().index();
        int j = bond.second().index();
        char order = static_cast<char>(bond.order());
        bond_info_[i * n + j] = order;
        bond_info_[j * n + i] = order;
    }
}


void Molecule::init_bond_distances() {
    const size_t n = atoms_->size();
    bond_distances_.resize(n * n);
    std::fill(bond_distances_.begin(), bond_distances_.end(), -1);

    for (size_t i = 0; i < n; i++) {
        auto q = std::queue<int>();
        q.push(static_cast<int>(i));
        bond_distances_[i * n + i] = 0;
        while (!q.empty()) {
            int p = q.front();
            q.pop();
            for (int neighbor: get_bonded(p)) {
                if (bond_distances_[i * n + neighbor] == -1) {
                    q.push(neighbor);
                    bond_distances_[i * n + neighbor] = 1 + bond_distances_[i * n + p];
                }
            }
        }
    }
}


int Molecule::bond_distance(const Atom &atom1, const Atom &atom2) const {
    const size_t n = atoms_->size();
    return bond_distances_[atom1.index() * n + atom2.index()];
}


std::vector<const Atom *> Molecule::k_bond_distance(const Atom &atom, int k) const {
    const size_t n = atoms_->size();
    std::vector<const Atom *> res;
    for (size_t i = 0; i < n; i++) {
        if (bond_distances_[atom.index() * n + i] == k) {
            res.push_back(&atoms_->operator[](k));
        }
    }
    return res;
}


int Molecule::total_charge() const {
    int sum = 0;
    for (const auto &atom: *atoms_) {
        sum += atom.formal_charge();
    }
    return sum;
}


std::vector<const Atom *> Molecule::get_close_atoms(const Atom &atom, double cutoff) const {
    std::vector<const Atom *> atoms;
    for (const auto &other_atom: *atoms_) {
        if (other_atom == atom) {
            atoms.insert(atoms.begin(), &atom);
            continue;
        }
        if (distance(atom, other_atom) < cutoff) {
            atoms.push_back(&other_atom);
        }
    }
    return atoms;
}
