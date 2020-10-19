//
// Created by krab1k on 31/10/18.
//

#include <vector>
#include <cmath>
#include <functional>
#include <Eigen/LU>

#include "smpqeq.h"
#include "../parameters.h"
#include "../geometry.h"

CHARGEFW2_METHOD(SMP_QEq)


Eigen::VectorXd SMP_QEq::EE_system(const std::vector<const Atom *> &atoms, double total_charge) const {

    size_t n = atoms.size();

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + 1);

    for (int iter = 0; iter < 5; iter++)
    {
        for (size_t i = 0; i < n; i++) {
            const auto &atom_i = *atoms[i];
            A(i, i) = 2 * (parameters_->atom()->parameter(atom::second)(atom_i) +
                                parameters_->atom()->parameter(atom::third)(atom_i) * b(i) +
                                parameters_->atom()->parameter(atom::fourth)(atom_i) * b(i) * b(i));
            b(i) = -parameters_->atom()->parameter(atom::first)(atom_i);
            for (size_t j = i + 1; j < n; j++) {
                const auto &atom_j = *atoms[j];
                auto gamma = 2 * std::sqrt(parameters_->atom()->parameter(atom::second)(atom_i) *
                                           parameters_->atom()->parameter(atom::second)(atom_j));
                auto expr = 1 / std::cbrt(1 / std::pow(gamma, 3) + std::pow(distance(atom_i, atom_j), 3));
                A(i, j) = expr;
                A(j, i) = expr;
            }
        }

        A.row(n) = Eigen::VectorXd::Constant(n + 1, 1);
        A.col(n) = Eigen::VectorXd::Constant(n + 1, 1);
        A(n, n) = 0;
        b(n) = total_charge;

        b = A.partialPivLu().solve(b);
    }

    return b.head(n);
}


std::vector<double> SMP_QEq::calculate_charges(const Molecule &molecule) const {
    auto f = [this](const std::vector<const Atom *> &atoms, double total_charge) -> Eigen::VectorXd {
        return EE_system(atoms, total_charge);
    };

    Eigen::VectorXd q = solve_EE(molecule, f);
    return std::vector<double>(q.data(), q.data() + q.size());
}
