#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class SQE : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "SQE",
        .internal_name = "sqe",
        .full_name = "Split-charge equilibration",
        .publication = "10.1063/1.2346671",
        .type = "3D",
        .priority = 90
    };

    enum atom{electronegativity, hardness, width};
    enum bond{kappa};

public:
    explicit SQE() : Method({}, {"electronegativity", "hardness", "width"}, {"kappa"}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
