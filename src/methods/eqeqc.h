//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"


class EQeqC : public EEMethod {
    enum common{alpha};
    enum atom{Dz};
    std::vector<double> solve_system(const std::vector<const Atom *> &atoms, double total_charge) const override;
public:
    explicit EQeqC() : EEMethod("EQeq+C", {"alpha"}, {"Dz"}, {}, {}) {}

    virtual ~EQeqC() = default;
};

extern "C" BOOST_SYMBOL_EXPORT EQeqC method;
EQeqC method;
