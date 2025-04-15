#pragma once

#include <Eigen/Core>
#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class EEM : public EEMethod {
    enum common {kappa};
    enum atom {A, B};

    [[nodiscard]] Eigen::VectorXd EE_system(const std::vector<const Atom *> &atoms, double total_charge) const;

public:
    explicit EEM() : EEMethod({"kappa"}, {"A", "B"}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& metadata() const override;
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
