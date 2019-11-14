//
// Created by krab1k on 31/10/18.
//

#include <vector>
#include <Eigen/LU>

#include "tsef.h"
#include "../parameters.h"

CHARGEFW2_METHOD(TSEF)


double K(int i);


double K(int i) {
    double vals[] = {0.556, 0.778, 1.000, 1.053, 1.087, 1.091};
    if (i > 6)
        return vals[5];
    else
        return vals[i - 1];
}


std::vector<double> TSEF::calculate_charges(const Molecule &molecule) const {

    size_t n = molecule.atoms().size();

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + 1);

    const double alpha = 14.4;

    for (size_t i = 0; i < n; i++) {
        const auto &atom_i = molecule.atoms()[i];
        A(i, i) = parameters_->atom()->parameter(atom::hardness)(atom_i);
        b(i) = - parameters_->atom()->parameter(atom::electronegativity)(atom_i);
        for (size_t j = i + 1; j < n; j++) {
            const auto &atom_j = molecule.atoms()[j];
            int bd = molecule.bond_distance(atom_i, atom_j);
            auto x = alpha * K(bd) / (0.84 * bd + 0.46);
            A(i, j) = x;
            A(j, i) = x;
        }
    }

    A.row(n) = Eigen::VectorXd::Constant(n + 1, 1);
    A.col(n) = Eigen::VectorXd::Constant(n + 1, 1);
    A(n, n) = 0;
    b(n) = molecule.total_charge();

    Eigen::VectorXd q = A.partialPivLu().solve(b).head(n);
    return std::vector<double>(q.data(), q.data() + q.size());
}
