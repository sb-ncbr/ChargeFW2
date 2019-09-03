//
// Created by krab1k on 13.11.18.
//

#pragma once

#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"


class GDAC : public Method {
    enum atom{A, B};
public:
    explicit GDAC() : Method("GDAC", {}, {"A", "B"}, {},
            {
                {"iters", {"iters", "Number of iterations", "int", "7", {}}}
            }) {}

    virtual ~GDAC() = default;

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};

extern "C" BOOST_SYMBOL_EXPORT GDAC method;
GDAC method;
