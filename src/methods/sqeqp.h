#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class SQEqp : public Method {
    enum atom{electronegativity, hardness, width, q0};
    enum bond{kappa};

public:
    explicit SQEqp() : Method({}, {"electronegativity", "hardness", "width", "q0"}, {"kappa"}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override;

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
