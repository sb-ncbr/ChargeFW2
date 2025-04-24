#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class Charge2 : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "Charge2",
        .internal_name = "charge2",
        .full_name = "Charge2",
        .publication = "10.1002/jcc.540030316",
        .type = "2D",
        .priority = 30
    };

    enum common{a1, a2, a3, b, c, alpha};
    enum atom{chi, P0, q0};

public:
    explicit Charge2() : Method({"a1", "a2", "a3", "b", "c", "alpha"}, {"chi", "P0", "q0"}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] std::vector<RequiredFeatures> get_requirements() const override {
        return {RequiredFeatures::BOND_DISTANCES};
    }
};
