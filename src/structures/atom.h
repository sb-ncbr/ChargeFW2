//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <array>
#include <string>
#include <tuple>
#include <fmt/format.h>

#include "../element.h"


class Molecule;

class Atom {
    size_t index_{};
    const Element *element_{};
    std::array<double, 3> pos_{};
    const Molecule *molecule_{};
    size_t type_{};
    int formal_charge_{};
    std::string atom_name_{};
    int residue_id_{};
    std::string residue_{};
    std::string chain_id_{};
    bool hetatm_{};

    friend class Molecule;

    friend class MoleculeSet;

public:
    Atom(size_t index, const Element *element, double x, double y, double z, std::string atom_name, int residue_id,
         std::string residue, std::string chain_id, bool hetatm);

    [[nodiscard]] size_t index() const { return index_; }

    [[nodiscard]] int formal_charge() const { return formal_charge_; }

    [[nodiscard]] const Element &element() const { return *element_; }

    [[nodiscard]] const Molecule *molecule() const { return molecule_; }

    [[nodiscard]] const std::array<double, 3> &pos() const { return pos_; }

    [[nodiscard]] size_t type() const { return type_; }

    [[nodiscard]] int residue_id() const { return residue_id_; }

    [[nodiscard]] std::string residue() const { return residue_; }

    [[nodiscard]] std::string chain_id() const { return chain_id_; }

    [[nodiscard]] std::string name() const { return atom_name_; }

    [[nodiscard]] bool hetatm() const { return hetatm_; }

    void _set_formal_charge(int charge) { formal_charge_ = charge; }

    bool inline operator==(const Atom &other) const {
        return this->index_ == other.index_ and this->molecule_ == other.molecule_;
    }
};


namespace fmt {
    template<>
    struct formatter<Atom> {
        template<typename ParseContext>
        constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }

        template<typename FormatContext>
        auto format(const Atom &a, FormatContext &ctx) {
            return format_to(ctx.begin(), "Atom {} Idx: {}", a.name(), a.index());
        }
    };
}
