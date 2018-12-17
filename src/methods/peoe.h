//
// Created by krab1k on 13.11.18.
//

#pragma once


#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"

class PEOE : public Method {
    enum common{dampH};
    enum atom{A, B, C};
public:
    explicit PEOE() : Method("PEOE", {"dampH"}, {"A", "B", "C"}, {},
            {
                {"iters", {"iters", "Number of iterations", "int", "7", {}}}
            }) {}

    std::vector<double> calculate_charges(const Molecule &molecule) override;
};

extern "C" BOOST_SYMBOL_EXPORT PEOE method;
PEOE method;