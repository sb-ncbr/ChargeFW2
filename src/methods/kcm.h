#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class KCM : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "KCM",
        .internal_name = "kcm",
        .full_name = "Kirchhoff Charge Model",
        .publication = "10.1002/jcc.20892",
        .type = "2D",
        .priority = 60
    };

    enum atom{electronegativity, hardness};
public:
    explicit KCM() : Method({}, {"electronegativity", "hardness"}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
