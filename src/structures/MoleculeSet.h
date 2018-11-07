//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <vector>
#include <memory>
#include <utility>
#include <string>

#include "Molecule.h"
#include "../Parameters.h"

class MoleculeSet {
    std::unique_ptr<std::vector<Molecule> > molecules_;
public:
    explicit MoleculeSet(std::unique_ptr<std::vector<Molecule> > molecules);

    void info() const;

    const std::vector<Molecule> &molecules() const { return *molecules_; }

    void classify_atoms(std::string classifier);
    void classify_atoms_from_parameters(const Parameters &parameters);
};
