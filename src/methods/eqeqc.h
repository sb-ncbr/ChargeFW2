//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <Eigen/Core>
#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"


class EQeqC : public EEMethod {
    enum common{alpha};
    enum atom{Dz};
    [[nodiscard]] Eigen::VectorXd solve_system(const std::vector<const Atom *> &atoms, double total_charge) const override;

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
public:
    explicit EQeqC() : EEMethod("EQeq+C", {"alpha"}, {"Dz"}, {}, {}) {}

    virtual ~EQeqC() = default;
};

extern "C" BOOST_SYMBOL_EXPORT EQeqC method;
EQeqC method;
