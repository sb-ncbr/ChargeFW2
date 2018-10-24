//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <vector>
#include <utility>

#include "Molecule.h"

class MoleculeSet {
    std::vector<Molecule> molecules;
public:
    explicit MoleculeSet(std::vector<Molecule> molecules) : molecules{std::move(molecules)} {};

    void info();
};


