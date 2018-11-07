#include <utility>

//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <Eigen/Dense>

#include "structures/Molecule.h"
#include "Parameters.h"

class Method {
protected:
    const std::vector<std::string> common_parameters_{};
    const std::vector<std::string> atom_parameters_{};
    const std::vector<std::string> bond_parameters_{};

    const Parameters *parameters_{nullptr};

public:
    Method(std::vector<std::string> common, std::vector<std::string> atom, std::vector<std::string> bond,
           const Parameters *parameters) :
            common_parameters_{std::move(common)},
            atom_parameters_{std::move(atom)},
            bond_parameters_{std::move(bond)} {
        set_parameters(parameters);
    }

    void set_parameters(const Parameters *parameters);

    virtual Eigen::VectorXd calculate_charges(const Molecule &molecule) = 0;
};

