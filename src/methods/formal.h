#pragma once

#include "../method.h"


class Formal : public Method {
public:
    explicit Formal() : Method("Formal charges (from file)", {}, {}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& get_metadata() const override;

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
