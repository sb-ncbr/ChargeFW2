#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class ABEEM : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "ABEEM",
        .internal_name = "abeem",
        .full_name = "Atom-Bond Electronegativity Equalization Method",
        .publication = "10.1021/jp9711048",
        .type = "3D",
        .priority = 190
    };

    enum common{k};
    enum atom{a, b, c};
    enum bond{A, B, C, D};

public:
    explicit ABEEM() : Method({"k"}, {"a", "b", "c"}, {"A", "B", "C", "D"}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
