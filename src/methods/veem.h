#pragma once

#include "../method.h"


class VEEM : public Method {
public:
    inline static const MethodMetadata METADATA = {
        .name = "VEEM",
        .internal_name = "veem",
        .full_name = "Valence Electrons Equalization Method",
        .publication = "10.1088/1674-0068/24/01/31-39",
        .type = "2D",
        .priority = 20
    };

    explicit VEEM() : Method({}, {}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_molecule(const Molecule &molecule) const override;
};
