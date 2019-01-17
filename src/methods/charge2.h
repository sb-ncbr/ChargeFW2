//
// Created by krab1k on 13.11.18.
//

#pragma once


#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"

class Charge2 : public Method {
    enum common{a1, a2, a3, b, c, alpha};
    enum atom{chi, P0, q0};
public:
    explicit Charge2() : Method("Charge2", {"a1", "a2", "a3", "b", "c", "alpha"}, {"chi", "P0", "q0"}, {}, {}) {}

    std::vector<double> calculate_charges(const Molecule &molecule) const override;
};

extern "C" BOOST_SYMBOL_EXPORT Charge2 method;
Charge2 method;