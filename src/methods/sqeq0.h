//
// Created by krab1k on 27.08.20.
//

#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class SQEq0 : public Method {
    enum atom{electronegativity, hardness, width};
    enum bond{kappa};

public:
    explicit SQEq0() : Method("SQEq0", {}, {"electronegativity", "hardness", "width"}, {"kappa"}, {}) {}

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
