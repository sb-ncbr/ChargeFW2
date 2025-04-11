#pragma once

#include <Eigen/Core>
#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class SFKEEM : public EEMethod {
    enum common{sigma};
    enum atom{A, B};
    [[nodiscard]] Eigen::VectorXd EE_system(const std::vector<const Atom *> &atoms, double total_charge) const;

public:
    explicit SFKEEM() : EEMethod("SFKEEM", {"sigma"}, {"A", "B"}, {}, {}) {}

    [[nodiscard]] const MethodMetadata& get_metadata() const override;
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
