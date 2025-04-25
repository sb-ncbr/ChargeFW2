#pragma once

#include <Eigen/Core>
#include <vector>

#include "../structures/molecule.h"
#include "../method.h"

class EQeq : public EEMethod {
    inline static const MethodMetadata METADATA = {
        .name = "EQeq",
        .internal_name = "eqeq",
        .full_name = "Extended Charge Equilibration Method",
        .publication = "10.1021/jz3008485",
        .type = "3D",
        .priority = 150
    };

    [[nodiscard]] Eigen::VectorXd EE_system(const std::vector<const Atom *> &atoms, double total_charge) const;

public:
    explicit EQeq() : EEMethod({}, {}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
