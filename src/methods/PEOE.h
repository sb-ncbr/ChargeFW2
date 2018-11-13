//
// Created by krab1k on 13.11.18.
//

#pragma once


#include <vector>
#include <boost/config.hpp>

#include "../structures/Molecule.h"
#include "../Method.h"

class PEOE : public Method {
    enum common{dampH};
    enum atom{A, B, C};
public:
    explicit PEOE() : Method({"dampH"}, {"A", "B", "C"}, {}) {}

    std::vector<double> calculate_charges(const Molecule &molecule) override;
};

extern "C" BOOST_SYMBOL_EXPORT PEOE method;
PEOE method;