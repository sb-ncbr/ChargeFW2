//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"

class EQeqC : public Method {
    enum common{alpha};
    enum atom{Dz};
public:
    explicit EQeqC() : Method("EQeq+C", {"alpha"}, {"Dz"}, {}, {}) {}

    std::vector<double> calculate_charges(const Molecule &molecule) override;
};

extern "C" BOOST_SYMBOL_EXPORT EQeqC method;
EQeqC method;