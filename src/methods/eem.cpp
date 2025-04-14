#include <vector>
#include <Eigen/LU>

#include "eem.h"
#include "../parameters.h"
#include "../geometry.h"

CHARGEFW2_METHOD(EEM)

namespace {
    const MethodMetadata EEM_METADATA = {
        .name = "EEM",
        .internal_name = "eem",
        .full_name = "Electronegativity Equalization Method",
        .publication = "10.1021/ja00275a013",
        .type = "3D",
        .priority = 200
    };
};

const MethodMetadata& EEM::metadata() const {
    return EEM_METADATA;
};

Eigen::VectorXd EEM::EE_system(const std::vector<const Atom *> &atoms, double total_charge) const {

    const auto n = static_cast<Eigen::Index>(atoms.size());

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + 1);

    for (Eigen::Index i = 0; i < n; i++) {
        const auto &atom_i = *atoms[i];
        A(i, i) = parameters_->atom()->parameter(atom::B)(atom_i);
        b(i) = -parameters_->atom()->parameter(atom::A)(atom_i);
        for (Eigen::Index j = i + 1; j < n; j++) {
            const auto &atom_j = *atoms[j];
            auto x = parameters_->common()->parameter(common::kappa) / distance(atom_i, atom_j);
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


std::vector<double> EEM::calculate_charges(const Molecule &molecule) const {
    auto f = [this](const std::vector<const Atom *> &atoms, double total_charge) -> Eigen::VectorXd {
        return EE_system(atoms, total_charge);
    };

    Eigen::VectorXd q = solve_EE(molecule, f);
    return {q.data(), q.data() + q.size()};
}
