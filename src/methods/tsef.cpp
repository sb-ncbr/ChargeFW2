//
// Created by krab1k on 31/10/18.
//

#include <vector>
#include <cmath>
#include <mkl_lapacke.h>
#include <mkl.h>

#include "tsef.h"
#include "../parameters.h"


#define IDX(i, j) (i * m + j)


double K(int i) {
    double vals[] = {0.556, 0.778, 1.000, 1.053, 1.087, 1.091};
    if (i > 6)
        return vals[5];
    else
        return vals[i - 1];
}


std::vector<double> TSEF::calculate_charges(const Molecule &molecule) {

    size_t n = molecule.atoms().size();
    size_t m = n + 1;

    auto *A = (double *) mkl_malloc(m * m * sizeof(double), 64);
    auto *b = (double *) mkl_malloc(m * sizeof(double), 64);
    auto *ipiv = (MKL_INT *) mkl_malloc(m * sizeof(MKL_INT), 64);

    const double alpha = 14.4;

    for (size_t i = 0; i < n; i++) {
        const auto &atom_i = molecule.atoms()[i];
        A[IDX(i, i)] = parameters_->atom()->parameter(atom::hardness)(atom_i);
        b[i] = - parameters_->atom()->parameter(atom::electronegativity)(atom_i);
        for (size_t j = i + 1; j < n; j++) {
            const auto &atom_j = molecule.atoms()[j];
            int bd = molecule.bond_distance(atom_i, atom_j);
            A[IDX(i, j)] = alpha * K(bd) / (0.84 * bd + 0.46);
        }
    }

    for(size_t i = 0; i < n; i++) {
        A[IDX(i, n)] = 1;
    }

    A[IDX(n, n)] = 0;
    b[n] = molecule.total_charge();

    MKL_INT info = LAPACKE_dsysv(LAPACK_ROW_MAJOR, 'U', m, 1, A, m, ipiv, b, 1);
    if(info) {
        throw std::runtime_error("Cannot solve linear system");
    }

    std::vector<double> results;
    results.assign(b, b + n);

    mkl_free(A);
    mkl_free(b);
    mkl_free(ipiv);

    return results;
}