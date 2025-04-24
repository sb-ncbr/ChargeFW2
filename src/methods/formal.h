#pragma once

#include "../method.h"


class Formal : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "Formal charges (from file)",
        .internal_name = "formal",
        .full_name = "Formal charges",
        .publication = std::nullopt,
        .type = "other",
        .priority = 10
    };
public:
    explicit Formal() : Method({}, {}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
