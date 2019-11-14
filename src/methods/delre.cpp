//
// Created by krab1k on 13.11.18.
//

#include <cmath>
#include <Eigen/LU>

#include "delre.h"
#include "../structures/molecule.h"
#include "../structures/bond.h"
#include "../parameters.h"

CHARGEFW2_METHOD(DelRe)


std::vector<double> DelRe::calculate_charges(const Molecule &molecule) const {

    const size_t n = molecule.atoms().size();
    const size_t m = molecule.bonds().size();

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n);

    for (size_t i = 0; i < n; i++) {
        auto &atom_i = molecule.atoms()[i];
        b(i) = -parameters_->atom()->parameter(atom::delta)(atom_i);
        A(i, i) = -1.0;
    }

    for (const auto &bond: molecule.bonds()) {
        size_t i = bond.first().index();
        size_t j = bond.second().index();
        A(i, j) = parameters_->bond()->parameter(bond::gammaA)(bond);
        A(j, i) = parameters_->bond()->parameter(bond::gammaB)(bond);
    }

    Eigen::VectorXd d = A.partialPivLu().solve(b);
    std::vector<double> q(n, 0);

    for (size_t k = 0; k < m; k++) {
        const auto &bond = molecule.bonds()[k];
        size_t i = bond.first().index();
        size_t j = bond.second().index();
        double dq = (d(i) - d(j)) / (2 * parameters_->bond()->parameter(bond::eps)(bond));
        q[i] -= dq;
        q[j] += dq;
    }

    return q;
}
