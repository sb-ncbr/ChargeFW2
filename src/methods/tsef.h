//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"


class TSEF : public Method {
    enum atom{electronegativity, hardness};
public:
    explicit TSEF() : Method("TSEF", {}, {"electronegativity", "hardness"}, {}, {}) {}

    virtual ~TSEF() = default;

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] std::vector<RequiredFeatures> get_requirements() const override {
        return {RequiredFeatures::BOND_DISTANCES};
    }
};

extern "C" BOOST_SYMBOL_EXPORT TSEF method;
TSEF method;
