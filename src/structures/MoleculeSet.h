//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <vector>
#include <memory>
#include <utility>

#include "Molecule.h"

class MoleculeSet {
    std::unique_ptr<std::vector<Molecule> > molecules_;
public:
    explicit MoleculeSet(std::unique_ptr<std::vector<Molecule> > molecules);

    void info() const;

    const std::vector<Molecule> &molecules() const { return *molecules_; }

    void classify_atoms(QString classifier);
};
