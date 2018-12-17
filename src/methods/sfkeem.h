//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"

class SFKEEM : public Method {
    enum common{sigma};
    enum atom{A, B};
public:
    explicit SFKEEM() : Method("SFKEEM", {"sigma"}, {"A", "B"}, {}, {}) {}

    std::vector<double> calculate_charges(const Molecule &molecule) override;
};

extern "C" BOOST_SYMBOL_EXPORT SFKEEM method;
SFKEEM method;