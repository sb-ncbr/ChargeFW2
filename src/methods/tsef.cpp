#include <vector>
#include <Eigen/LU>

#include "tsef.h"
#include "../parameters.h"

CHARGEFW2_METHOD(TSEF)

namespace {
    const MethodMetadata TSEF_METADATA = {
        .name = "TSEF",
        .internal_name = "tsef",
        .full_name = "Topologically Symmetrical Energy Function",
        .publication = "10.1080/10629360701844142",
        .type = "2D",
        .priority = 40
    };
};

const MethodMetadata& TSEF::metadata() const {
    return TSEF_METADATA;
};

double K(int i);


double K(int i) {
    double vals[] = {0.556, 0.778, 1.000, 1.053, 1.087, 1.091};
    if (i > 6)
        return vals[5];
    else
        return vals[i - 1];
}


std::vector<double> TSEF::calculate_charges(const Molecule &molecule) const {

    const auto n = static_cast<Eigen::Index>(molecule.atoms().size());

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + 1);

    const double alpha = 14.4;

    for (Eigen::Index i = 0; i < n; i++) {
        const auto &atom_i = molecule.atoms()[i];
        A(i, i) = parameters_->atom()->parameter(atom::hardness)(atom_i);
        b(i) = - parameters_->atom()->parameter(atom::electronegativity)(atom_i);
        for (Eigen::Index j = i + 1; j < n; j++) {
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
    return {q.data(), q.data() + q.size()};
}
