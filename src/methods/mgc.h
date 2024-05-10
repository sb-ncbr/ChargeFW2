#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class MGC : public Method {
public:
    explicit MGC() : Method("MGC", {}, {}, {}, {}) {}

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
