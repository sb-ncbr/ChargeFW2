#pragma once

#include <Eigen/Core>
#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class SMP_QEq : public EEMethod {
    inline static const MethodMetadata METADATA = {
        .name = "SMP/QEq",
        .internal_name = "smpqeq",
        .full_name = "Self-Consistent Charge Equilibration Method",
        .publication = "10.1021/jp8063273",
        .type = "3D",
        .priority = 160
    };

    enum atom{first, second, third, fourth};

    [[nodiscard]] Eigen::VectorXd EE_system(const std::vector<const Atom *> &atoms, double total_charge) const;

public:
    explicit SMP_QEq() : EEMethod({}, {"first", "second", "third", "fourth"}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
