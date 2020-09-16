//
// Created by krab1k on 27.08.20.
//

#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class SQEqp : public Method {
    enum atom{electronegativity, hardness, width, q0};
    enum bond{kappa};

public:
    explicit SQEqp() : Method("SQE+qp", {}, {"electronegativity", "hardness", "width", "q0"}, {"kappa"}, {}) {}

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
