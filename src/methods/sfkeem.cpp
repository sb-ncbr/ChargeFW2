#include <vector>
#include <cmath>
#include <functional>
#include <Eigen/LU>

#include "sfkeem.h"
#include "../parameters.h"
#include "../geometry.h"

CHARGEFW2_METHOD(SFKEEM)

namespace {
    const MethodMetadata SFKEEM_METADATA = {
        .internal_name = "sfkeem",
        .full_name = "Selfconsistent Functional Kernel Equalized Electronegativity Method",
        .publication = "10.1021/ci050505e",
        .type = "3D",
        .priority = 180,
        .has_parameters = true
    };
};

const MethodMetadata& SFKEEM::get_metadata() const {
    return SFKEEM_METADATA;
};

Eigen::VectorXd SFKEEM::EE_system(const std::vector<const Atom *> &atoms, double total_charge) const {

    const auto n = static_cast<Eigen::Index>(atoms.size());

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + 1);

    for (Eigen::Index i = 0; i < n; i++) {
        const auto &atom_i = *atoms[i];
        A(i, i) = 2 * parameters_->atom()->parameter(atom::B)(atom_i);
        b(i) = -parameters_->atom()->parameter(atom::A)(atom_i);
        for (Eigen::Index j = i + 1; j < n; j++) {
            const auto &atom_j = *atoms[j];
            auto x = 2 * sqrt(parameters_->atom()->parameter(atom::B)(atom_i) *
                                    parameters_->atom()->parameter(atom::B)(atom_j)) /
                           cosh(parameters_->common()->parameter(common::sigma) * distance(atom_i, atom_j));
            A(i, j) = x;
            A(j, i) = x;
        }
    }

    A.row(n) = Eigen::VectorXd::Constant(n + 1, 1);
    A.col(n) = Eigen::VectorXd::Constant(n + 1, 1);
    A(n, n) = 0;
    b(n) = total_charge;

    return A.partialPivLu().solve(b).head(n);
}


std::vector<double> SFKEEM::calculate_charges(const Molecule &molecule) const {
    auto f = [this](const std::vector<const Atom *> &atoms, double total_charge) -> Eigen::VectorXd {
        return EE_system(atoms, total_charge);
    };

    Eigen::VectorXd q = solve_EE(molecule, f);
    return {q.data(), q.data() + q.size()};
}
