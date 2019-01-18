//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"

class SMP_QEq : public EEMethod {
    enum atom{first, second, third, fourth};
    std::vector<double> solve_system(const std::vector<const Atom *> &atoms, double total_charge) const override;
public:
    explicit SMP_QEq() : EEMethod("SMP/QEq", {}, {"first", "second", "third", "fourth"}, {}, {}) {}
};

extern "C" BOOST_SYMBOL_EXPORT SMP_QEq method;
SMP_QEq method;