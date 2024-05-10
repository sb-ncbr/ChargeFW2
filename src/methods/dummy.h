#pragma once

#include "../method.h"


class Dummy : public Method {
public:
    explicit Dummy() : Method("Dummy", {}, {}, {}, {}) {}

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
