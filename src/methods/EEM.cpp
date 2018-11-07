//
// Created by krab1k on 31/10/18.
//

#include "EEM.h"

#include <Eigen/Core>
#include <vector>
#include <cmath>


double distance(const Atom &atom1, const Atom &atom2) {
    double dx = atom1.pos()[0] - atom2.pos()[0];
    double dy = atom1.pos()[1] - atom2.pos()[1];
    double dz = atom1.pos()[2] - atom2.pos()[2];

    return std::sqrt(dx * dx + dy * dy + dz * dz);
}


Eigen::VectorXd EEM::calculate_charges(const Molecule &molecule) {

    long n = molecule.atoms().size();
    Eigen::MatrixXd matrix(n + 1, n + 1);
    Eigen::VectorXd vec(n + 1);

    for (int i = 0; i < n; i++) {
        const auto &atom_i = molecule.atoms()[i];
        matrix(i, i) = parameters_->atom()->parameter(atom::B)(atom_i);
        vec(i) = - parameters_->atom()->parameter(atom::A)(atom_i);
        for (int j = i + 1; j < n; j++) {
            const auto &atom_j = molecule.atoms()[j];
            matrix(i, j) = parameters_->common()->parameter(common::kappa) / distance(atom_i, atom_j);
            matrix(j, i) = matrix(i, j);
        }
    }

    matrix.row(n).setConstant(1);
    matrix.col(n).setConstant(1);
    matrix(n, n) = 0;
    vec(n) = 0;

    return matrix.partialPivLu().solve(vec).head(n);
}