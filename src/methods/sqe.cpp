#include <vector>
#include <cmath>
#include <Eigen/LU>

#include "sqe.h"
#include "../parameters.h"
#include "../geometry.h"

CHARGEFW2_METHOD(SQE)

namespace {
    const MethodMetadata SQE_METADATA = {
        .internal_name = "sqe",
        .full_name = "Split-charge equilibration",
        .publication = "10.1063/1.2346671",
        .type = "3D",
        .priority = 90,
        .has_parameters = true
    };
};

const MethodMetadata& SQE::get_metadata() const {
    return SQE_METADATA;
};

std::vector<double> SQE::calculate_charges(const Molecule &molecule) const {

    const auto n = static_cast<Eigen::Index>(molecule.atoms().size());
    const auto m = static_cast<Eigen::Index>(molecule.bonds().size());

    Eigen::VectorXd hardness = Eigen::VectorXd::Zero(n);

    for (Eigen::Index i = 0; i < n; i++) {
        const auto &atom = molecule.atoms()[i];
        hardness(i) = parameters_->atom()->parameter(atom::hardness)(atom);
    }

    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(m, n);
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(molecule.bonds().size()); i++) {
        const auto &bond = molecule.bonds()[i];
        auto i1 = static_cast<Eigen::Index>(bond.first().index());
        auto i2 = static_cast<Eigen::Index>(bond.second().index());
        T(i, i1) = 1;
        T(i, i2) = -1;
    }

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n);

    /* Setup EEM part */
    for (Eigen::Index i = 0; i < n; i++) {
        const auto &atom_i = molecule.atoms()[i];
        A(i, i) = hardness(i);
        b(i) = -parameters_->atom()->parameter(atom::electronegativity)(atom_i);
        for (Eigen::Index j = i + 1; j < n; j++) {
            const auto &atom_j = molecule.atoms()[j];
            auto d = distance(atom_i, atom_j);
            auto wi = parameters_->atom()->parameter(atom::width)(atom_i);
            auto wj = parameters_->atom()->parameter(atom::width)(atom_j);
            auto d0 = sqrt(2 * wi * wi + 2 * wj * wj);
            auto x = erf(d / d0) / d;
            A(i, j) = x;
            A(j, i) = x;
        }
    }

    Eigen::MatrixXd split_A = T * A * T.transpose();
    Eigen::VectorXd split_b = T * b;

    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(molecule.bonds().size()); i++) {
        const auto &bond = molecule.bonds()[i];
        split_A(i, i) += parameters_->bond()->parameter(bond::kappa)(bond);
    }

    Eigen::VectorXd split_q = split_A.partialPivLu().solve(split_b);
    Eigen::VectorXd q = T.transpose() * split_q;

    return {q.data(), q.data() + n};
}
