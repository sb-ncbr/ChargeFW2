//
// Created by krab1k on 23/10/18.
//

#include <queue>
#include <utility>
#include <tuple>
#include <map>
#include <string>
#include <vector>
#include <nanoflann.hpp>

#include "atom.h"
#include "bond.h"
#include "molecule.h"


Molecule::Molecule(std::string name, std::unique_ptr<std::vector<Atom> > atoms,
                   std::unique_ptr<std::vector<Bond> > bonds) {
    name_ = std::move(name);
    atoms_ = std::move(atoms);
    bonds_ = std::move(bonds);

    /* Calculate max bond orders */
    const size_t n = atoms_->size();
    max_hbo_.resize(n, 0);

    for (const auto &bond: *bonds_) {
        size_t i1 = bond.first().index_;
        size_t i2 = bond.second().index_;
        int bo = bond.order();
        max_hbo_[i1] = std::max(max_hbo_[i1], bo);
        max_hbo_[i2] = std::max(max_hbo_[i2], bo);
    }

    /* Calculate bonded elements */
    std::vector<std::vector<std::string>> neighbours;
    neighbours.resize(atoms_->size());
    for (const auto &bond: *bonds_) {
        const auto &atom1 = bond.first();
        const auto &atom2 = bond.second();
        neighbours[atom1.index()].emplace_back(atom2.element().symbol());
        neighbours[atom2.index()].emplace_back(atom1.element().symbol());
    }

    for (size_t i = 0; i < atoms_->size(); i++) {
        auto &bonded = neighbours[i];
        std::sort(bonded.begin(), bonded.end());

        std::string atom_type;
        for (const auto &symbol: bonded) {
            atom_type += symbol;
        }
        neighbour_elements_.push_back(atom_type);
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


std::vector<size_t> Molecule::get_bonded(size_t atom_idx) const {
    const size_t n = atoms_->size();
    std::vector<size_t> res;

    for (size_t j = 0; j < n; j++) {
        if (bond_info_[atom_idx * n + j]) {
            res.push_back(static_cast<size_t>(j));
        }
    }
    return res;
}


void Molecule::init_bond_info() {
    const size_t n = atoms_->size();
    bond_info_.resize(n * n);

    for (const auto &bond: *bonds_) {
        size_t i = bond.first().index();
        size_t j = bond.second().index();
        auto order = static_cast<char>(bond.order());
        bond_info_[i * n + j] = order;
        bond_info_[j * n + i] = order;
    }
}


void Molecule::init_bond_distances() {
    const size_t n = atoms_->size();
    bond_distances_.resize(n * n);
    std::fill(bond_distances_.begin(), bond_distances_.end(), -1);

    for (size_t i = 0; i < n; i++) {
        auto q = std::queue<size_t>();
        q.push(i);
        bond_distances_[i * n + i] = 0;
        while (!q.empty()) {
            size_t p = q.front();
            q.pop();
            for (size_t neighbor: get_bonded(p)) {
                if (bond_distances_[i * n + neighbor] == -1) {
                    q.push(neighbor);
                    bond_distances_[i * n + neighbor] = 1 + bond_distances_[i * n + p];
                }
            }
        }
    }
}


void Molecule::init_distance_tree() {
    adaptor_ = std::make_unique<AtomKDTreeAdaptor>(this);
    index_ = std::make_unique<kdtree_t>(3, *adaptor_, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    index_->buildIndex();
}


int Molecule::bond_distance(const Atom &atom1, const Atom &atom2) const {
    const size_t n = atoms_->size();
    return bond_distances_[atom1.index() * n + atom2.index()];
}


std::vector<const Atom *> Molecule::k_bond_distance(const Atom &atom, size_t k) const {
    const size_t n = atoms_->size();
    std::vector<const Atom *> res;
    for (size_t i = 0; i < n; i++) {
        if (bond_distances_[atom.index() * n + i] == static_cast<int>(k)) {
            res.push_back(&atoms_->operator[](i));
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
    std::vector<const Atom *> close_atoms;

    std::vector<std::pair<uint32_t , double>> results;
    nanoflann::SearchParams params;

    auto matches_count = index_->radiusSearch(atom.pos().data(), cutoff * cutoff, results, params);
    for (size_t i = 0; i < matches_count; i++) {
        if (results[i].first == atom.index()) {
            close_atoms.insert(close_atoms.begin(), &atom);
            continue;
        }
        close_atoms.push_back(&atoms()[results[i].first]);
    }
    return close_atoms;
}


const Bond *Molecule::get_bond(const Atom &atom1, const Atom &atom2) const {
    for (const auto &bond: *bonds_) {
        if ((atom1 == bond.first() and atom2 == bond.second()) or (atom1 == bond.second() and atom2 == bond.first()))
            return &bond;
    }
    return nullptr;
}
