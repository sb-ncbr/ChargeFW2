#pragma once

#include "../method.h"


class Dummy : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "Dummy",
        .internal_name = "dummy",
        .full_name = "Dummy Method",
        .publication = std::nullopt,
        .type = "other",
        .priority = 0
    };
public:
    explicit Dummy() : Method({}, {}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
