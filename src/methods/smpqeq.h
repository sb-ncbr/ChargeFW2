//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <Eigen/Core>
#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"


class SMP_QEq : public EEMethod {
    enum atom{first, second, third, fourth};
    [[nodiscard]] Eigen::VectorXd solve_system(const std::vector<const Atom *> &atoms, double total_charge) const override;

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
public:
    explicit SMP_QEq() : EEMethod("SMP/QEq", {}, {"first", "second", "third", "fourth"}, {}, {}) {}

    virtual ~SMP_QEq() = default;
};

extern "C" BOOST_SYMBOL_EXPORT SMP_QEq method;
SMP_QEq method;
