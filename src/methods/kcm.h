#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class KCM : public Method {
    enum atom{electronegativity, hardness};
public:
    explicit KCM() : Method({}, {"electronegativity", "hardness"}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override;
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
