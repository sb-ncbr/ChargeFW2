#include <vector>
#include <Eigen/LU>

#include "denr.h"
#include "../parameters.h"

CHARGEFW2_METHOD(DENR)

namespace {
    const MethodMetadata DENR_METADATA = {
        .name = "DENR",
        .internal_name = "denr",
        .full_name = "Dynamical Electronegativity Relaxation",
        .publication = "10.1080/10629360701844142",
        .type = "2D",
        .priority = 50
    };
};

const MethodMetadata& DENR::metadata() const {
    return DENR_METADATA;
};

std::vector<double> DENR::calculate_charges(const Molecule &molecule) const {

    const auto n = static_cast<Eigen::Index>(molecule.atoms().size());

    Eigen::MatrixXd eta = Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd chi = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd q = Eigen::VectorXd::Zero(n);

    for (Eigen::Index i = 0; i < n; i++) {
        auto &atom_i = molecule.atoms()[i];
        chi(i) = parameters_->atom()->parameter(atom::electronegativity)(atom_i);
        eta(i, i) = parameters_->atom()->parameter(atom::hardness)(atom_i);
    }

    for (const auto &bond: molecule.bonds()) {
        auto i1 = static_cast<Eigen::Index>(bond.first().index());
        auto i2 = static_cast<Eigen::Index>(bond.second().index());
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

    return {q.data(), q.data() + q.size()};
}
