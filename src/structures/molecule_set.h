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
#include "../method.h"
#include "../parameters.h"


enum class AtomClassifier {
    PLAIN,
    HBO
};


enum class BondClassifier {
    PLAIN,
    BO
};


class MoleculeSet {
    std::vector<std::tuple<std::string, std::string, std::string>> atom_types_{};
    std::vector<std::tuple<std::string, std::string, std::string, std::string>> bond_types_{};
    std::unique_ptr<std::vector<Molecule> > molecules_{nullptr};

    size_t classify_atoms_from_parameters(const Parameters &parameters, bool remove_unclassified = true,
                                          bool permissive_types = false);

    size_t classify_bonds_from_parameters(const Parameters &parameters, bool remove_unclassified = true,
                                          bool permissive_types = false);

    [[nodiscard]] std::vector<int> get_max_bond_orders(const Molecule &molecule) const;

    void set_atom_type(Atom &atom, const std::tuple<std::string, std::string, std::string> &tuple);

    void set_bond_type(Bond &bond, const std::tuple<std::string, std::string, std::string, std::string> &tuple);

public:
    explicit MoleculeSet(std::unique_ptr<std::vector<Molecule> > molecules);

    void info() const;

    [[nodiscard]] const std::vector<Molecule> &molecules() const { return *molecules_; }

    void classify_atoms(AtomClassifier cls);

    void classify_bonds(BondClassifier cls);

    size_t classify_set_from_parameters(const Parameters &parameters, bool remove_unclassified = true,
                                        bool permissive_types = false);

    [[nodiscard]] std::vector<std::tuple<std::string, std::string, std::string>>
    atom_types() const { return atom_types_; }

    [[nodiscard]] std::vector<std::tuple<std::string, std::string, std::string, std::string>>
    bond_types() const { return bond_types_; }

    void fulfill_requirements(const std::vector<RequiredFeatures> &features);

    [[nodiscard]] bool has_proteins() const;
};
