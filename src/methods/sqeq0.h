#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class SQEq0 : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "SQEq0",
        .internal_name = "sqeq0",
        .full_name = "Split-charge equilibration with initial formal charges",
        .publication = "10.1021/ct200512e",
        .type = "3D",
        .priority = 80
    };

    enum atom{electronegativity, hardness, width};
    enum bond{kappa};

public:
    explicit SQEq0() : Method({}, {"electronegativity", "hardness", "width"}, {"kappa"}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
