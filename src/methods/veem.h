#pragma once

#include "../method.h"


class VEEM : public Method {
public:
    explicit VEEM() : Method("VEEM", {}, {}, {}, {}) {}

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_molecule(const Molecule &molecule) const override;
};
