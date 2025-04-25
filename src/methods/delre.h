#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class DelRe : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "DelRe",
        .internal_name = "delre",
        .full_name = "Method of Del Re",
        .publication = "10.1039/JR9580004031",
        .type = "2D",
        .priority = 130
    };

    enum atom{delta};
    enum bond{eps, gammaA, gammaB};

public:
    explicit DelRe() : Method({}, {"delta"}, {"eps", "gammaA", "gammaB"}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
