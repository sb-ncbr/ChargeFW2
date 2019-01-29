//
// Created by krab1k on 31/10/18.
//

#include <vector>
#include <cmath>
#include <mkl_lapacke.h>
#include <mkl.h>

#include "sfkeem.h"
#include "../parameters.h"
#include "../geometry.h"

#define IDX(i, j) ((i) * m + (j))


std::vector<double> SFKEEM::solve_system(const std::vector<const Atom *> &atoms, double total_charge) const {

    size_t n = atoms.size();
    size_t m = n + 1;

    auto *A = static_cast<double *>(mkl_malloc(m * m * sizeof(double), 64));
    auto *b = static_cast<double *>(mkl_malloc(m * sizeof(double), 64));
    auto *ipiv = static_cast<int *>(mkl_malloc(m * sizeof(int), 64));

    for (size_t i = 0; i < n; i++) {
        const auto &atom_i = *atoms[i];
        A[IDX(i, i)] = 2 * parameters_->atom()->parameter(atom::B)(atom_i);
        b[i] = -parameters_->atom()->parameter(atom::A)(atom_i);
        for (size_t j = i + 1; j < n; j++) {
            const auto &atom_j = *atoms[j];
            A[IDX(i, j)] = 2 * sqrt(parameters_->atom()->parameter(atom::B)(atom_i) *
                                    parameters_->atom()->parameter(atom::B)(atom_j)) /
                           cosh(parameters_->common()->parameter(common::sigma) * distance(atom_i, atom_j));
        }
    }

    for (size_t i = 0; i < n; i++) {
        A[IDX(i, n)] = 1;
    }

    A[IDX(n, n)] = 0;
    b[n] = total_charge;

    int info = LAPACKE_dsysv(LAPACK_ROW_MAJOR, 'U', m, 1, A, m, ipiv, b, 1);
    if (info) {
        throw std::runtime_error("Cannot solve linear system");
    }

    std::vector<double> results;
    results.assign(b, b + n);

    mkl_free(A);
    mkl_free(b);
    mkl_free(ipiv);

    return results;
}
