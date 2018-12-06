//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <vector>
#include <memory>
#include <utility>
#include <tuple>
#include <string>

#include "molecule.h"
#include "../classifier.h"
#include "../parameters.h"

class Classifier;

class MoleculeSet {
    std::vector<std::tuple<std::string, std::string, std::string>> atom_types_{};
    std::unique_ptr<std::vector<Molecule> > molecules_{nullptr};
public:
    explicit MoleculeSet(std::unique_ptr<std::vector<Molecule> > molecules);

    void info() const;

    const std::vector<Molecule> &molecules() const { return *molecules_; }

    void classify_atoms(const Classifier &cls);
    void classify_atoms_from_parameters(const Parameters &parameters);
    int get_unclassified_molecules_count(const Parameters &parameters) const;
};
