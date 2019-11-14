//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <Eigen/Core>
#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class SMP_QEq : public EEMethod {
    enum atom{first, second, third, fourth};

    [[nodiscard]] Eigen::VectorXd EE_system(const std::vector<const Atom *> &atoms, double total_charge) const;

public:
    explicit SMP_QEq() : EEMethod("SMP/QEq", {}, {"first", "second", "third", "fourth"}, {}, {}) {}

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
