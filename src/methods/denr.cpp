//
// Created by krab1k on 31/10/18.
//

#include <vector>
#include <cmath>
#include <mkl_lapacke.h>
#include <mkl.h>

#include "denr.h"
#include "../structures/molecule.h"
#include "../parameters.h"


std::vector<double> DENR::calculate_charges(const Molecule &molecule) const {

    size_t n = molecule.atoms().size();

    auto *L = static_cast<double *>(mkl_calloc(n * n, sizeof(double), 64));
    auto *chi = static_cast<double *>(mkl_calloc(n, sizeof(double), 64));
    auto *q = static_cast<double *>(mkl_calloc(n, sizeof(double), 64));
    auto *ipiv = static_cast<int *>(mkl_calloc(n, sizeof(int), 64));
    auto *tmp = static_cast<double *>(mkl_calloc(n, sizeof(double), 64));
    auto *I = static_cast<double *>(mkl_calloc(n * n, sizeof(double), 64));
    auto *eta = static_cast<double *>(mkl_calloc(n * n, sizeof(double), 64));

    for (size_t i = 0; i < n; i++) {
        auto &atom_i = molecule.atoms()[i];
        chi[i] = parameters_->atom()->parameter(atom::electronegativity)(atom_i);
        eta[i * n + i] = parameters_->atom()->parameter(atom::hardness)(atom_i);
        I[i * n + i] = 1.0;
        L[i * n + i] = molecule.degree(atom_i);
        for (size_t j = i + 1; j < n; j++) {
            auto &atom_j = molecule.atoms()[j];
            if (molecule.bonded(atom_i, atom_j)) {
                L[i * n + j] = -1;
                L[j * n + i] = -1;
            }
        }
    }

    double step = parameters_->common()->parameter(common::step);

    auto n_int = static_cast<int>(n);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_int, n_int, n_int, step, L, n_int, eta, n_int, 1.0, I, n_int);
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n_int, n_int, I, n_int, ipiv);

    cblas_dgemv(CblasRowMajor, CblasNoTrans, n_int, n_int, step, L, n_int, chi, 1, 0, tmp, 1);

    for (int i = 0; i < parameters_->common()->parameter(common::iterations); i++) {
        cblas_daxpby(n_int, -1.0, tmp, 1, 1.0, q, 1);
        LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', n_int, 1, I, n_int, ipiv, q, 1);
    }

    std::vector<double> res(n, 0);
    for (size_t i = 0; i < n; i++) {
        res[i] = q[i];
    }

    mkl_free(L);
    mkl_free(chi);
    mkl_free(ipiv);
    mkl_free(q);
    mkl_free(eta);
    mkl_free(tmp);
    mkl_free(I);

    return res;
}
