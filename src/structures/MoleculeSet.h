//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <QVector>
#include <utility>

#include "Molecule.h"

class MoleculeSet {
    QVector<Molecule> molecules;
public:
    explicit MoleculeSet(const QVector<Molecule> &molecules) : molecules{std::move(molecules)} {};

    void info();
};


