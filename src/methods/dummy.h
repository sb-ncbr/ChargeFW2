#pragma once

#include "../method.h"


class Dummy : public Method {
public:
    explicit Dummy() : Method("Dummy", {}, {}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& get_metadata() const override;

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
