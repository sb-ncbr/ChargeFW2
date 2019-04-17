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

    virtual ~DENR() = default;

    std::vector<double> calculate_charges(const Molecule &molecule) const override;

    bool is_suitable_for_large_molecule() const override { return false; }
};

extern "C" BOOST_SYMBOL_EXPORT DENR method;
DENR method;
