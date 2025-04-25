#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class DENR : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "DENR",
        .internal_name = "denr",
        .full_name = "Dynamical Electronegativity Relaxation",
        .publication = "10.1080/10629360701844142",
        .type = "2D",
        .priority = 50
    };

    enum common{step, iterations};
    enum atom{electronegativity, hardness};
public:
    explicit DENR() : Method({"step", "iterations"}, {"electronegativity", "hardness"}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
