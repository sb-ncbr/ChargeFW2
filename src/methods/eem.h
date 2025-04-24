#pragma once

#include <Eigen/Core>
#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class EEM : public EEMethod {
    inline static const MethodMetadata METADATA = {
        .name = "EEM",
        .internal_name = "eem",
        .full_name = "Electronegativity Equalization Method",
        .publication = "10.1021/ja00275a013",
        .type = "3D",
        .priority = 200
    };

    enum common {kappa};
    enum atom {A, B};

    [[nodiscard]] Eigen::VectorXd EE_system(const std::vector<const Atom *> &atoms, double total_charge) const;

public:
    explicit EEM() : EEMethod({"kappa"}, {"A", "B"}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
