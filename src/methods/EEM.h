//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <Eigen/Core>

#include "../structures/Molecule.h"
#include "../Method.h"

class EEM : public Method {
public:
    explicit EEM(const Parameters *parameters) : Method({"kappa"}, {"A", "B"}, {}, parameters) {}

    Eigen::VectorXd calculate_charges(const Molecule &molecule) override;
};
