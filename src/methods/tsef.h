//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"

class TSEF : public Method {
    enum atom{electronegativity, hardness};
public:
    explicit TSEF() : Method("TSEF", {}, {"electronegativity", "hardness"}, {}, {}) {}

    std::vector<double> calculate_charges(const Molecule &molecule) const override;
};

extern "C" BOOST_SYMBOL_EXPORT TSEF method;
TSEF method;