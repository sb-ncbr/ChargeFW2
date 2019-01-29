//
// Created by krab1k on 8.11.18.
//

#pragma once

#include <boost/config.hpp>

#include "../method.h"


class VEEM : public Method {
public:
    explicit VEEM() : Method("VEEM", {}, {}, {}, {}) {}

    std::vector<double> calculate_charges(const Molecule &molecule) const override;
};

extern "C" BOOST_SYMBOL_EXPORT VEEM method;
VEEM method;
