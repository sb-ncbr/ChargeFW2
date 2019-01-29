//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"


class SFKEEM : public EEMethod {
    enum common{sigma};
    enum atom{A, B};
    std::vector<double> solve_system(const std::vector<const Atom *> &atoms, double total_charge) const override;
public:
    explicit SFKEEM() : EEMethod("SFKEEM", {"sigma"}, {"A", "B"}, {}, {}) {}


};

extern "C" BOOST_SYMBOL_EXPORT SFKEEM method;
SFKEEM method;
