#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class DENR : public Method {
    enum common{step, iterations};
    enum atom{electronegativity, hardness};
public:
    explicit DENR() : Method("DENR", {"step", "iterations"}, {"electronegativity", "hardness"}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& get_metadata() const override;

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
