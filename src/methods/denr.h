//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"

class DENR : public Method {
    enum common{step, iterations};
    enum atom{electronegativity, hardness};
public:
    explicit DENR() : Method("DENR", {"step", "iterations"}, {"electronegativity", "hardness"}, {}, {}) {}

    std::vector<double> calculate_charges(const Molecule &molecule) const override;

    std::vector<RequiredFeatures> get_requirements() const override {
        return {RequiredFeatures::BOND_INFO};
    }
};

extern "C" BOOST_SYMBOL_EXPORT DENR method;
DENR method;