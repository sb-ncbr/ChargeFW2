//
// Created by krab1k on 13.11.18.
//

#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class MPEOE : public Method {
    enum common{Hplus};
    enum atom{A, B};
    enum bond{f};
public:
    explicit MPEOE() : Method("MPEOE", {"Hplus"}, {"A", "B"}, {"f"},
            {
                {"iters", {"iters", "Number of iterations", "int", "7", {}}}
            }
    ) {}

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
