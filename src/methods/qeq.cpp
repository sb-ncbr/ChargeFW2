//
// Created by krab1k on 02/01/19.
//

#include <functional>
#include <string>
#include <cmath>
#include <Eigen/LU>

#include "qeq.h"
#include "../structures/atom.h"
#include "../geometry.h"
#include "../parameters.h"

CHARGEFW2_METHOD(QEq)


double QEq::overlap_term(const Atom &atom_i, const Atom &atom_j, const std::string &type) const {
    auto Ji = parameters_->atom()->parameter(atom::hardness)(atom_i);
    auto Jj = parameters_->atom()->parameter(atom::hardness)(atom_j);
    auto Rij = distance(atom_i, atom_j);
    if (type == "Nishimoto-Mataga") {
        return 1 / (Rij + 2 / (Ji + Jj));
    } else if (type == "Nishimoto-Mataga-Weiss") {
        const double f = 1.2;
        return f / (Rij + (2 * f) / (Ji + Jj));
    } else if (type == "Ohno") {
        return 1 / std::sqrt(Rij * Rij + std::pow(2 / (Ji + Jj), 2));
    } else if (type == "Ohno-Klopman") {
        return 1 / std::sqrt(Rij * Rij + std::pow(1 / (2 * Ji) + 1 / (2 * Jj), 2));
    } else if (type == "DasGupta-Huzinaga") {
        const double k = 0.4;
        return 1 / (Rij + 1 / (Ji / 2 * exp(k * Rij) + Jj / 2 * exp(k * Rij)));
    } else /* (type == "Louwen-Vogt") */ {
        const double gamma = (Ji + Jj) / 2;
        return 1 / std::cbrt(1 / std::pow(gamma, 3) + std::pow(Rij, 3));
    }
}


Eigen::VectorXd QEq::EE_system(const std::vector<const Atom *> &atoms, double total_charge) const {

    size_t n = atoms.size();

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + 1);

    const auto type = get_option_value<std::string>("overlap_term");

    for (size_t i = 0; i < n; i++) {
        const auto &atom_i = *atoms[i];
        A(i, i) = parameters_->atom()->parameter(atom::hardness)(atom_i);
        b(i) = - parameters_->atom()->parameter(atom::electronegativity)(atom_i);
        for (size_t j = i + 1; j < n; j++) {
            const auto &atom_j = *atoms[j];
            auto x = overlap_term(atom_i, atom_j, type);
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


std::vector<double> QEq::calculate_charges(const Molecule &molecule) const {
    auto f = [this](const std::vector<const Atom *> &atoms, double total_charge) -> Eigen::VectorXd {
        return EE_system(atoms, total_charge);
    };

    Eigen::VectorXd q = solve_EE(molecule, f);
    return std::vector<double>(q.data(), q.data() + q.size());
}
