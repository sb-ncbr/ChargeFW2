//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <QString>
#include <iostream>
#include <tuple>
#include <Eigen/Dense>

#include "../Element.h"

class Molecule;

class Atom {
    int index_{};
    const Element *element_{};
    Eigen::Vector3d pos_{};
    const Molecule *molecule_{};
    std::tuple<QString, QString, QString> atom_type_{};

    friend class Molecule;

    friend class MoleculeSet;

public:
    Atom(int index, const Element *element, double x, double y, double z);

    int index() const { return index_; }

    const Element &element() const { return *element_; }

    const Molecule *molecule() const { return molecule_; }

    const Eigen::Vector3d &pos() const { return pos_; }

    const std::tuple<QString, QString, QString> atom_type() const { return atom_type_; }

    friend std::ostream &operator<<(std::ostream &str, const Atom &atom);

    bool operator==(const Atom &other) const;
};
