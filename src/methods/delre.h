#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class DelRe : public Method {
    enum atom{delta};
    enum bond{eps, gammaA, gammaB};

public:
    explicit DelRe() : Method({}, {"delta"}, {"eps", "gammaA", "gammaB"}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override;
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
