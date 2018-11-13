//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <iostream>

#include "Atom.h"

class Bond {
    const Atom *first_{};
    const Atom *second_{};
    int order_{};
    const Molecule *molecule_{};
public:
    Bond(const Atom *atom1, const Atom *atom2, int order) : first_{atom1}, second_{atom2}, order_{order} {}

    bool hasAtom(const Atom &atom) const { return atom == *first_ or atom == *second_; }

    const Atom &first() const { return *first_; }

    const Atom &second() const { return *second_; }

    int order() const { return order_; }

    friend std::ostream &operator<<(std::ostream &str, const Bond &bond);

    friend class MoleculeSet;
};
