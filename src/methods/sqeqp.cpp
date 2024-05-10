//
// Created by krab1k on 27.08.20.
//

#include <vector>
#include <cmath>
#include <Eigen/LU>

#include "sqeqp.h"
#include "../parameters.h"
#include "../geometry.h"

CHARGEFW2_METHOD(SQEqp)


std::vector<double> SQEqp::calculate_charges(const Molecule &molecule) const {

    const auto n = static_cast<Eigen::Index>(molecule.atoms().size());
    const auto m = static_cast<Eigen::Index>(molecule.bonds().size());

    Eigen::VectorXd q0 = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd hardness = Eigen::VectorXd::Zero(n);

    for (Eigen::Index i = 0; i < n; i++) {
        const auto &atom = molecule.atoms()[i];
        q0(i) = parameters_->atom()->parameter(atom::q0)(atom);
        hardness(i) = parameters_->atom()->parameter(atom::hardness)(atom);
    }

    q0 = q0.array() - (q0.sum() - molecule.total_charge()) / static_cast<double>(n);

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

    b = b - A * q0;
    b = b + hardness.cwiseProduct(q0);

    Eigen::MatrixXd split_A = T * A * T.transpose();
    Eigen::VectorXd split_b = T * b;

    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(molecule.bonds().size()); i++) {
        const auto &bond = molecule.bonds()[i];
        split_A(i, i) += parameters_->bond()->parameter(bond::kappa)(bond);
    }

    Eigen::VectorXd split_q = split_A.partialPivLu().solve(split_b);
    Eigen::VectorXd q = T.transpose() * split_q + q0;

    return {q.data(), q.data() + n};
}
