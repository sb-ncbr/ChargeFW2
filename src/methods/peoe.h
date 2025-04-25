#pragma once


#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class PEOE : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "PEOE",
        .internal_name = "peoe",
        .full_name = "Partial Equalization of Atomic Electronegativity",
        .publication = "10.1016/0040-4020(80)80168-2",
        .type = "2D",
        .priority = 120
    };

    enum common{dampH};
    enum atom{A, B, C};
public:
    explicit PEOE() : Method({"dampH"}, {"A", "B", "C"}, {},
            {
                {"iters", {"iters", "Number of iterations", "int", "7", {}}}
            }) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
