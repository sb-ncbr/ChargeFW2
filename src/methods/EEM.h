//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <boost/config.hpp>

#include "../structures/Molecule.h"
#include "../Method.h"

class EEM : public Method {
    enum common{kappa};
    enum atom{A, B};
public:
    explicit EEM() : Method({"kappa"}, {"A", "B"}, {}) {}

    std::vector<double> calculate_charges(const Molecule &molecule) override;
};

extern "C" BOOST_SYMBOL_EXPORT EEM method;
EEM method;