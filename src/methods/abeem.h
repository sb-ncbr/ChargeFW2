//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class ABEEM : public Method {
    enum common{k};
    enum atom{a, b, c};
    enum bond{A, B, C, D};

public:
    explicit ABEEM() : Method("ABEEM", {"k"}, {"a", "b", "c"}, {"A", "B", "C", "D"}, {}) {}

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
