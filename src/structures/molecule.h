//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <fmt/format.h>
#include <utility>
#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <memory>
#include <nanoflann.hpp>

#include "atom.h"
#include "bond.h"


class AtomKDTreeAdaptor;


typedef nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, AtomKDTreeAdaptor>, AtomKDTreeAdaptor, 3> kdtree_t;


class Molecule {
    std::string name_;
    std::unique_ptr<std::vector<Atom> > atoms_;
    std::unique_ptr<std::vector<Bond> > bonds_;
    std::vector<char> bond_info_{};
    std::vector<int> bond_distances_{};
    std::unique_ptr<kdtree_t> index_{nullptr};
    std::unique_ptr<AtomKDTreeAdaptor> adaptor_{nullptr};

    std::vector<size_t> get_bonded(size_t atom_idx) const;

    void init_bond_info();

    void init_bond_distances();

    void init_distance_tree();

public:
    const std::vector<Atom> &atoms() const { return *atoms_; }

    const std::vector<Bond> &bonds() const { return *bonds_; }

    const std::string &name() const { return name_; }

    bool bonded(const Atom &atom1, const Atom &atom2) const;

    int bond_order(const Atom &atom1, const Atom &atom2) const;

    int degree(const Atom &atom) const;

    std::vector<const Atom *> get_close_atoms(const Atom &atom, double cutoff) const;

    Molecule() = default;

    Molecule(std::string name, std::unique_ptr<std::vector<Atom> > atoms, std::unique_ptr<std::vector<Bond> > bonds,
             const std::map<size_t, int> &charges);

    int bond_distance(const Atom &atom1, const Atom &atom2) const;

    std::vector<const Atom *> k_bond_distance(const Atom &atom, size_t k) const;

    int total_charge() const;

    friend class MoleculeSet;
};


namespace fmt {
    template<>
    struct formatter<Molecule> {
        template<typename ParseContext>
        constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }

        template<typename FormatContext>
        auto format(const Molecule &m, FormatContext &ctx) {
            return format_to(ctx.begin(), "Molecule {} Atoms: {} Bonds {}\n", m.name(), m.atoms().size(),
                             m.bonds().size());
        }
    };
}


class AtomKDTreeAdaptor {
    const Molecule *molecule_;
public:
    AtomKDTreeAdaptor(const Molecule *molecule) : molecule_{molecule} {};

    inline size_t kdtree_get_point_count() const {
        return molecule_->atoms().size();
    }

    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return molecule_->atoms()[idx].pos()[dim];
    }

    template<class BBOX>
    bool kdtree_get_bbox(BBOX &) const { return false; }
};
