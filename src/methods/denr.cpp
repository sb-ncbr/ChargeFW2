//
// Created by krab1k on 31/10/18.
//

#include <vector>
#include <cmath>
#include <Eigen/LU>

#include "denr.h"
#include "../structures/molecule.h"
#include "../parameters.h"

CHARGEFW2_METHOD(DENR)


std::vector<double> DENR::calculate_charges(const Molecule &molecule) const {

    size_t n = molecule.atoms().size();

    Eigen::MatrixXd eta = Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd chi = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd q = Eigen::VectorXd::Zero(n);

    for (size_t i = 0; i < n; i++) {
        auto &atom_i = molecule.atoms()[i];
        chi(i) = parameters_->atom()->parameter(atom::electronegativity)(atom_i);
        eta(i, i) = parameters_->atom()->parameter(atom::hardness)(atom_i);
    }

    for (const auto &bond: molecule.bonds()) {
        auto i1 = bond.first().index();
        auto i2 = bond.second().index();
        L(i1, i1) += 1;
        L(i2, i2) += 1;
        L(i1, i2) -= 1;
        L(i2, i1) -= 1;
    }

    double step = parameters_->common()->parameter(common::step);

    Eigen::PartialPivLU<Eigen::MatrixXd> x = (Eigen::MatrixXd::Identity(n, n) + step * L * eta).partialPivLu();
    Eigen::VectorXd tmp = step * L * chi;
    for (int i = 0; i < parameters_->common()->parameter(common::iterations); i++) {
        q = x.solve(q - tmp);
    }

    return std::vector<double>(q.data(), q.data() + q.size());
}
