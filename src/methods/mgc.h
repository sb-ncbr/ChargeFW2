#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class MGC : public Method {
    inline static const MethodMetadata METADATA = {
        .name = "MGC",
        .internal_name = "mgc",
        .full_name = "Molecular Graph Charge",
        .publication = "10.1002/poc.378",
        .type = "2D",
        .priority = 70
    };

public:
    explicit MGC() : Method({}, {}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    [[nodiscard]] bool is_suitable_for_large_molecule() const override { return false; }
};
