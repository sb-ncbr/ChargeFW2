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
    const std::vector<QString> common_parameters_{};
    const std::vector<QString> atom_parameters_{};
    const std::vector<QString> bond_parameters_{};

    const Parameters *parameters_{nullptr};

public:
    Method(std::vector<QString> common, std::vector<QString> atom, std::vector<QString> bond,
           const Parameters *parameters) :
            common_parameters_{std::move(common)},
            atom_parameters_{std::move(atom)},
            bond_parameters_{std::move(bond)} {
        set_parameters(parameters);
    }

    void set_parameters(const Parameters *parameters);

    virtual Eigen::VectorXd calculate_charges(const Molecule &molecule) = 0;
};

