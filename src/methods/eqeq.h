//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <Eigen/Core>
#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"


class EQeq : public EEMethod {
    [[nodiscard]] Eigen::VectorXd solve_system(const std::vector<const Atom *> &atoms, double total_charge) const override;

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
public:
    explicit EQeq() : EEMethod("EQeq", {}, {}, {}, {}) {}

    virtual ~EQeq() = default;
};

extern "C" BOOST_SYMBOL_EXPORT EQeq method;
EQeq method;
