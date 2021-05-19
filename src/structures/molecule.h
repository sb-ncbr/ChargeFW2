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
    std::vector<int> max_hbo_{};
    std::vector<std::string> neighbour_elements_{};
    std::vector<char> bond_info_{};
    std::vector<int> bond_distances_{};
    std::unique_ptr<kdtree_t> index_{nullptr};
    std::unique_ptr<AtomKDTreeAdaptor> adaptor_{nullptr};

    [[nodiscard]] std::vector<size_t> get_bonded(size_t atom_idx) const;

    void init_bond_info();

    void init_bond_distances();

    void init_distance_tree();

public:
    [[nodiscard]] const std::vector<Atom> &atoms() const { return *atoms_; }

    [[nodiscard]] const std::vector<Bond> &bonds() const { return *bonds_; }

    [[nodiscard]] const std::string &name() const { return name_; }

    [[nodiscard]] bool bonded(const Atom &atom1, const Atom &atom2) const;

    [[nodiscard]] const Bond *get_bond(const Atom &atom1, const Atom &atom2) const;

    [[nodiscard]] int bond_order(const Atom &atom1, const Atom &atom2) const;

    [[nodiscard]] int degree(const Atom &atom) const;

    [[nodiscard]] const std::vector<int> &get_max_bond_orders() const { return max_hbo_; }

    [[nodiscard]] const std::vector<std::string> &get_bonded_elements() const { return neighbour_elements_; }

    [[nodiscard]] std::vector<const Atom *> get_close_atoms(const Atom &atom, double cutoff) const;

    Molecule() = default;

    Molecule(std::string name, std::unique_ptr<std::vector<Atom> > atoms, std::unique_ptr<std::vector<Bond> > bonds);

    [[nodiscard]] int bond_distance(const Atom &atom1, const Atom &atom2) const;

    [[nodiscard]] std::vector<const Atom *> k_bond_distance(const Atom &atom, size_t k) const;

    [[nodiscard]] int total_charge() const;

    [[nodiscard]] bool is_protein() const { return not(*atoms_)[0].chain_id().empty(); }

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
    explicit AtomKDTreeAdaptor(const Molecule *molecule) : molecule_{molecule} {}

    [[nodiscard]] inline size_t kdtree_get_point_count() const {
        return molecule_->atoms().size();
    }

    [[nodiscard]] inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return molecule_->atoms()[idx].pos()[dim];
    }

    template<class BBOX>
    bool kdtree_get_bbox(BBOX &) const { return false; }
};
