#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class TSEF : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "TSEF",
        .internal_name = "tsef",
        .full_name = "Topologically Symmetrical Energy Function",
        .publication = "10.1080/10629360701844142",
        .type = "2D",
        .priority = 40
    };

    enum atom{electronegativity, hardness};

public:
    explicit TSEF() : Method({}, {"electronegativity", "hardness"}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] std::vector<RequiredFeatures> get_requirements() const override {
        return {RequiredFeatures::BOND_DISTANCES};
    }
};
