#pragma once

#include <Eigen/Core>
#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class EQeqC : public EEMethod {
    inline static const MethodMetadata METADATA = {
        .name = "EQeq+C",
        .internal_name = "eqeqc",
        .full_name = "Bond-Order-Corrected Extended Charge Equilibration Method",
        .publication = "10.1021/acs.jctc.5b00037",
        .type = "3D",
        .priority = 140
    };

    enum common{alpha};
    enum atom{Dz};

    [[nodiscard]] Eigen::VectorXd EE_system(const std::vector<const Atom *> &atoms, double total_charge) const;

public:
    explicit EQeqC() : EEMethod({"alpha"}, {"Dz"}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
