#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class Charge2 : public Method {
    enum common{a1, a2, a3, b, c, alpha};
    enum atom{chi, P0, q0};

public:
    explicit Charge2() : Method("Charge2", {"a1", "a2", "a3", "b", "c", "alpha"}, {"chi", "P0", "q0"}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& get_metadata() const override;
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] std::vector<RequiredFeatures> get_requirements() const override {
        return {RequiredFeatures::BOND_DISTANCES};
    }
};
