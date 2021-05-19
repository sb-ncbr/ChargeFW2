//
// Created by krab1k on 31/10/18.
//

#include <functional>
#include <vector>
#include <cmath>
#include <Eigen/LU>

#include "eqeq.h"
#include "../parameters.h"
#include "../geometry.h"

CHARGEFW2_METHOD(EQeq)


Eigen::VectorXd EQeq::EE_system(const std::vector<const Atom *> &atoms, double total_charge) const {

    size_t n = atoms.size();

    const double lambda = 1.2;
    const double k = 14.4;
    double H_electron_affinity = -2.0; // Exception for hydrogen mentioned in the article

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + 1);
    Eigen::VectorXd J = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd X = Eigen::VectorXd::Zero(n);

    for (size_t i = 0; i < n; i++) {
        const auto &atom_i = *atoms[i];
        if (atom_i.element().symbol() == "H") {
            X(i) = (atom_i.element().ionization_potential() + H_electron_affinity) / 2;
            J(i) = atom_i.element().ionization_potential() - H_electron_affinity;
        } else {
            X(i) = (atom_i.element().ionization_potential() + atom_i.element().electron_affinity()) / 2;
            J(i) = atom_i.element().ionization_potential() - atom_i.element().electron_affinity();
        }
    }

    for (size_t i = 0; i < n; i++) {
        const auto &atom_i = *atoms[i];
        A(i, i) = J(i);
        b(i) = -X(i);
        for (size_t j = i + 1; j < n; j++) {
            const auto &atom_j = *atoms[j];
            double a = std::sqrt(J(i) * J(j)) / k;
            double Rij = distance(atom_i, atom_j);
            double overlap = std::exp(-a * a * Rij * Rij) * (2 * a - a * a * Rij - 1 / Rij);
            auto x = lambda * k / 2 * (1 / Rij + overlap);
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


std::vector<double> EQeq::calculate_charges(const Molecule &molecule) const {
    auto f = [this](const std::vector<const Atom *> &atoms, double total_charge) -> Eigen::VectorXd {
        return EE_system(atoms, total_charge);
    };

    Eigen::VectorXd q = solve_EE(molecule, f);

    return std::vector<double>(q.data(), q.data() + q.size());
}
