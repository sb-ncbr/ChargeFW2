//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"

class SMP_QEq : public Method {
    enum atom{first, second, third, fourth};
public:
    explicit SMP_QEq() : Method("SMP/QEq", {}, {"first", "second", "third", "fourth"}, {}, {}) {}

    std::vector<double> calculate_charges(const Molecule &molecule) const override;
};

extern "C" BOOST_SYMBOL_EXPORT SMP_QEq method;
SMP_QEq method;