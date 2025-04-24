#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class SQEqp : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "SQE+qp",
        .internal_name = "sqeqp",
        .full_name = "Split-charge equilibration with parametrized initial charges",
        .publication = "10.1186/s13321-021-00528-w",
        .type = "3D",
        .priority = 210
    };

    enum atom{electronegativity, hardness, width, q0};
    enum bond{kappa};

public:
    explicit SQEqp() : Method({}, {"electronegativity", "hardness", "width", "q0"}, {"kappa"}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
