//
// Created by krab1k on 13.11.18.
//

#include <cmath>
#include <mkl_lapacke.h>
#include <mkl.h>

#include "delre.h"
#include "../structures/molecule.h"
#include "../structures/bond.h"
#include "../parameters.h"


std::vector<double> DelRe::calculate_charges(const Molecule &molecule) {

    const size_t n = molecule.atoms().size();
    const size_t m = molecule.bonds().size();
    std::vector<double> q(n, 0);

    auto *A = (double *) mkl_calloc(n * n, sizeof(double), 64);
    auto *b = (double *) mkl_calloc(n, sizeof(double), 64);
    auto *ipiv = (lapack_int *) mkl_calloc(n, sizeof(lapack_int), 64);

    for (size_t i = 0; i < n; i++) {
        auto &atom_i = molecule.atoms()[i];
        b[i] = -parameters_->atom()->parameter(atom::delta)(atom_i);
        A[i * n + i] = -1.0;
    }

    for (const auto &bond: molecule.bonds()) {
        int i = bond.first().index();
        int j = bond.second().index();
        A[i * n + j] = parameters_->bond()->parameter(bond::gammaA)(bond);
        A[j * n + i] = parameters_->bond()->parameter(bond::gammaB)(bond);
    }

    LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, 1, A, n, ipiv, b, 1);

    for (size_t k = 0; k < m; k++) {
        const auto &bond = molecule.bonds()[k];
        int i = bond.first().index();
        int j = bond.second().index();
        double dq = (b[i] - b[j]) / (2 * parameters_->bond()->parameter(bond::eps)(bond));
        q[i] -= dq;
        q[j] += dq;
    }

    mkl_free(A);
    mkl_free(b);
    mkl_free(ipiv);

    return q;
}