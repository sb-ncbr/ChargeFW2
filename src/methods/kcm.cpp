//
// Created by krab1k on 31/10/18.
//

#include <vector>
#include <cmath>
#include <Eigen/LU>

#include "kcm.h"
#include "../structures/molecule.h"
#include "../parameters.h"

CHARGEFW2_METHOD(KCM)


std::vector<double> KCM::calculate_charges(const Molecule &molecule) const {

    const auto n = static_cast<Eigen::Index>(molecule.atoms().size());
    const auto m = static_cast<Eigen::Index>(molecule.bonds().size());

    Eigen::MatrixXd W = Eigen::MatrixXd::Zero(m, m);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m, n);
    Eigen::VectorXd chi0 = Eigen::VectorXd::Zero(n);

    /* Compute
     *  q = (B.T @ W @ B + I)^-1 @ chi0 - chi0
     */

    for (Eigen::Index i = 0; i < n; i++) {
        chi0(i) = parameters_->atom()->parameter(atom::electronegativity)(molecule.atoms()[i]);
    }

    for (Eigen::Index i = 0; i < m; i++) {
        auto &bond = molecule.bonds()[i];
        auto &first = bond.first();
        auto &second = bond.second();

        W(i, i) = 1 / (parameters_->atom()->parameter(atom::hardness)(first) +
                            parameters_->atom()->parameter(atom::hardness)(second));

        B(i, static_cast<Eigen::Index>(first.index())) = 1;
        B(i, static_cast<Eigen::Index>(second.index())) = -1;
    }

    Eigen::VectorXd q = (B.transpose() * W * B + Eigen::MatrixXd::Identity(n, n)).partialPivLu().solve(chi0) - chi0;
    return {q.data(), q.data() + q.size()};
}
