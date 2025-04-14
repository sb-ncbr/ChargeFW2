#pragma once

#include "../method.h"


class Dummy : public Method {
public:
    explicit Dummy() : Method({}, {}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override;

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
