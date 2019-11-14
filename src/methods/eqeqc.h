//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <Eigen/Core>
#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class EQeqC : public EEMethod {
    enum common{alpha};
    enum atom{Dz};

    [[nodiscard]] Eigen::VectorXd EE_system(const std::vector<const Atom *> &atoms, double total_charge) const;

public:
    explicit EQeqC() : EEMethod("EQeq+C", {"alpha"}, {"Dz"}, {}, {}) {}

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
