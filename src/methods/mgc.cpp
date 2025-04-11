#include <vector>
#include <cmath>
#include <Eigen/LU>

#include "mgc.h"
#include "../structures/molecule.h"

CHARGEFW2_METHOD(MGC)

namespace {
    const MethodMetadata MGC_METADATA = {
        .internal_name = "mgc",
        .full_name = "Molecular Graph Charge",
        .publication = "10.1002/poc.378",
        .type = "2D",
        .priority = 70,
        .has_parameters = false
    };
};

const MethodMetadata& MGC::get_metadata() const {
    return MGC_METADATA;
};

std::vector<double> MGC::calculate_charges(const Molecule &molecule) const {

    const auto n = static_cast<Eigen::Index>(molecule.atoms().size());

    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd X0 = Eigen::VectorXd::Zero(n);

    double log_sum = 0;

    for (const auto &atom: molecule.atoms()) {
        auto i = static_cast<Eigen::Index>(atom.index());
        S(i, i) = 1;
        X0(i) = atom.element().electronegativity();
        log_sum += log(X0(i));
    }

    for (const auto &bond: molecule.bonds()) {
        auto i1 = static_cast<Eigen::Index>(bond.first().index());
        auto i2 = static_cast<Eigen::Index>(bond.second().index());
        auto order = bond.order();
        S(i1, i1) += order;
        S(i2, i2) += order;
        S(i1, i2) -= order;
        S(i2, i1) -= order;
    }

    Eigen::VectorXd chi = S.partialPivLu().solve(X0);
    for (Eigen::Index i = 0; i < n; i++) {
        chi(i) -= molecule.atoms()[i].element().electronegativity();
    }
    chi /= exp(log_sum / static_cast<double>(n));

    return {chi.data(), chi.data() + chi.size()};
}
