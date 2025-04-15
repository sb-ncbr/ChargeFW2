#pragma once

#include "../method.h"


class Formal : public Method {
public:
    explicit Formal() : Method({}, {}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override;

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
