//
// Created by krab1k on 13.11.18.
//

#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class DelRe : public Method {
    enum atom{delta};
    enum bond{eps, gammaA, gammaB};

public:
    explicit DelRe() : Method("DelRe", {}, {"delta"}, {"eps", "gammaA", "gammaB"}, {}) {}

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
