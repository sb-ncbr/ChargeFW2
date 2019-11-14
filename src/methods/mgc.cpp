//
// Created by krab1k on 31/10/18.
//

#include <vector>
#include <cmath>
#include <Eigen/LU>

#include "mgc.h"
#include "../structures/molecule.h"

CHARGEFW2_METHOD(MGC)


std::vector<double> MGC::calculate_charges(const Molecule &molecule) const {

    size_t n = molecule.atoms().size();

    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd X0 = Eigen::VectorXd::Zero(n);

    double log_sum = 0;

    for (const auto &atom: molecule.atoms()) {
        auto i = atom.index();
        S(i, i) = 1;
        X0(i) = atom.element().electronegativity();
        log_sum += log(X0(i));
    }

    for (const auto &bond: molecule.bonds()) {
        auto i1 = bond.first().index();
        auto i2 = bond.second().index();
        auto order = bond.order();
        S(i1, i1) += order;
        S(i2, i2) += order;
        S(i1, i2) -= order;
        S(i2, i1) -= order;
    }

    Eigen::VectorXd chi = S.partialPivLu().solve(X0);
    for (size_t i = 0; i < n; i++) {
        chi(i) -= molecule.atoms()[i].element().electronegativity();
    }
    chi /= exp(log_sum / n);

    return std::vector<double>(chi.data(), chi.data() + chi.size());
}
