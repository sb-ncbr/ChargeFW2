//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <array>
#include <string>
#include <iostream>
#include <tuple>

#include "../Element.h"

class Molecule;

class Atom {
    int index_{};
    const Element *element_{};
    std::array<double, 3> pos_{};
    const Molecule *molecule_{};
    std::tuple<std::string, std::string, std::string> atom_type_{};

    friend class Molecule;

    friend class MoleculeSet;

public:
    Atom(int index, const Element *element, double x, double y, double z);

    int index() const { return index_; }

    const Element &element() const { return *element_; }

    const Molecule *molecule() const { return molecule_; }

    const std::array<double, 3> &pos() const { return pos_; }

    const std::tuple<std::string, std::string, std::string> atom_type() const { return atom_type_; }

    friend std::ostream &operator<<(std::ostream &str, const Atom &atom);

    bool inline operator==(const Atom &other) const {
        return this->index_ == other.index_ and this->molecule_ == other.molecule_;
    };
};
