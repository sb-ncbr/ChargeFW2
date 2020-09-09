//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <array>

#include "atom.h"

class Bond {
    const Atom *first_{};
    const Atom *second_{};
    int order_{};
    size_t type_{};
    const Molecule *molecule_{};
public:
    Bond(const Atom *atom1, const Atom *atom2, int order) noexcept : first_{atom1}, second_{atom2}, order_{order} {}

    [[nodiscard]] bool hasAtom(const Atom &atom) const { return atom == *first_ or atom == *second_; }

    [[nodiscard]] const Atom &first() const { return *first_; }

    [[nodiscard]] const Atom &second() const { return *second_; }

    [[nodiscard]] int order() const { return order_; }

    [[nodiscard]] size_t type() const { return type_; }

    [[nodiscard]] std::array<double, 3> get_center(bool weighted=false) const;

    friend class MoleculeSet;
};


namespace fmt {
    template<>
    struct formatter<Bond> {
        template<typename ParseContext>
        constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }

        template<typename FormatContext>
        auto format(const Bond &b, FormatContext &ctx) {
            return format_to(ctx.begin(), "Bond ({}, {}): {}\n", b.first(), b.second(), b.order());
        }
    };
}
