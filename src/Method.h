#include <utility>

//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>

#include "structures/Molecule.h"
#include "Parameters.h"

class Method {
protected:
    const std::vector<std::string> common_parameters_{};
    const std::vector<std::string> atom_parameters_{};
    const std::vector<std::string> bond_parameters_{};

    const Parameters *parameters_{nullptr};

public:
    Method(std::vector<std::string> common, std::vector<std::string> atom, std::vector<std::string> bond) :
            common_parameters_{std::move(common)},
            atom_parameters_{std::move(atom)},
            bond_parameters_{std::move(bond)} {}

    void set_parameters(const Parameters *parameters);

    bool has_parameters() {
        return (common_parameters_.size() + atom_parameters_.size() + bond_parameters_.size()) != 0;
    }

    virtual std::vector<double> calculate_charges(const Molecule &molecule) = 0;
};

