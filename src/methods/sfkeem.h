#pragma once

#include <Eigen/Core>
#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class SFKEEM : public EEMethod {
    inline static const MethodMetadata METADATA = {
        .name = "SFKEEM",
        .internal_name = "sfkeem",
        .full_name = "Selfconsistent Functional Kernel Equalized Electronegativity Method",
        .publication = "10.1021/ci050505e",
        .type = "3D",
        .priority = 180
    };

    enum common{sigma};
    enum atom{A, B};
    [[nodiscard]] Eigen::VectorXd EE_system(const std::vector<const Atom *> &atoms, double total_charge) const;

public:
    explicit SFKEEM() : EEMethod({"sigma"}, {"A", "B"}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
