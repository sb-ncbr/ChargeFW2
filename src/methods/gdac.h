#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class GDAC : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "GDAC",
        .internal_name = "gdac",
        .full_name = "Geometry-Dependent Net Atomic Charges",
        .publication = "10.1021/jp0023213",
        .type = "3D",
        .priority = 100
    };

    enum atom{A, B};
public:
    explicit GDAC() : Method({}, {"A", "B"}, {},
            {
                {"iters", {"iters", "Number of iterations", "int", "7", {}}}
            }) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
