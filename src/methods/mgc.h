//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"


class MGC : public Method {
public:
    explicit MGC() : Method("MGC", {}, {}, {}, {}) {}

    virtual ~MGC() = default;

    std::vector<double> calculate_charges(const Molecule &molecule) const override;

    std::vector<RequiredFeatures> get_requirements() const override {
        return {RequiredFeatures::BOND_INFO};
    }

    bool is_suitable_for_large_molecule() const override { return false; }
};

extern "C" BOOST_SYMBOL_EXPORT MGC method;
MGC method;
