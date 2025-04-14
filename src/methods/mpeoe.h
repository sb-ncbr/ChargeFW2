#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class MPEOE : public Method {
    enum common{Hplus};
    enum atom{A, B};
    enum bond{f};
public:
    explicit MPEOE() : Method({"Hplus"}, {"A", "B"}, {"f"},
            {
                {"iters", {"iters", "Number of iterations", "int", "7", {}}}
            }
    ) {}

    [[nodiscard]] const MethodMetadata& metadata() const override;
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
