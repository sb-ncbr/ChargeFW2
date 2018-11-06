//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <array>
#include <QString>
#include <iostream>
#include <tuple>

#include "../Element.h"

class Molecule;

class Atom {
    int index_{};
    const Element *element_{};
    std::array<double, 3> pos_{};
    const Molecule *molecule_{};
    std::tuple<QString, QString, QString> atom_type_{};

public:
    Atom(int index, const Element *element, double x, double y, double z);

    int index() const { return index_; }

    const Element &element() const { return *element_; }

    const Molecule *molecule() const { return molecule_; }

    const std::tuple<QString, QString, QString> atom_type() const { return atom_type_; }

    friend std::ostream &operator<<(std::ostream &str, const Atom &atom);

    bool operator==(const Atom &other) const;

    friend class Molecule;

    friend class MoleculeSet;
};