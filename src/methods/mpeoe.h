#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class MPEOE : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "MPEOE",
        .internal_name = "mpeoe",
        .full_name = "Modified Partial Equalization of Atomic Electronegativity",
        .publication = "10.1021/j100374a066",
        .type = "2D",
        .priority = 110
    };

    enum common{Hplus};
    enum atom{A, B};
    enum bond{f};
public:
    explicit MPEOE() : Method({"Hplus"}, {"A", "B"}, {"f"},
            {
                {"iters", {"iters", "Number of iterations", "int", "7", {}}}
            }
    ) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
