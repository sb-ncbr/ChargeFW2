//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class TSEF : public Method {
    enum atom{electronegativity, hardness};

public:
    explicit TSEF() : Method("TSEF", {}, {"electronegativity", "hardness"}, {}, {}) {}

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] std::vector<RequiredFeatures> get_requirements() const override {
        return {RequiredFeatures::BOND_DISTANCES};
    }
};
