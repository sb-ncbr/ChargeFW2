//
// Created by krab1k on 31/10/18.
//

#include <vector>
#include <cmath>
#include <mkl_lapacke.h>
#include <mkl.h>

#include "kcm.h"
#include "../structures/molecule.h"
#include "../parameters.h"


std::vector<double> KCM::calculate_charges(const Molecule &molecule) const {

    size_t n = molecule.atoms().size();
    size_t m = molecule.bonds().size();

    auto *W = static_cast<double *>(mkl_calloc(m * m, sizeof(double), 64));
    auto *B = static_cast<double *>(mkl_calloc(m * n, sizeof(double), 64));
    auto *chi0 = static_cast<double *>(mkl_calloc(n, sizeof(double), 64));
    auto *tmp = static_cast<double *>(mkl_calloc(m * n, sizeof(double), 64));
    auto *res = static_cast<double *>(mkl_calloc(n * n, sizeof(double), 64));
    auto *ipiv = static_cast<int *>(mkl_malloc(n * sizeof(int), 64));

    /* Compute
     *  q = (B.T @ W @ B + I)^-1 @ chi0 - chi0
     */

    for (size_t i = 0; i < m; i++) {
        auto &bond = molecule.bonds()[i];
        auto &first = bond.first();
        auto &second = bond.second();

        W[m * i + i] = 1 / (parameters_->atom()->parameter(atom::hardness)(first) +
                            parameters_->atom()->parameter(atom::hardness)(second));

        B[i * n + first.index()] = 1;
        B[i * n + second.index()] = -1;
    }

    auto n_int = static_cast<int>(n);
    auto m_int = static_cast<int>(m);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m_int, n_int, m_int, 1.0, W, m_int, B, n_int, 0.0, tmp, n_int);

    std::vector<double> q(n, 0);
    for(size_t i = 0; i < n; i++) {
        res[i * n + i] = 1.0;
        double chi = parameters_->atom()->parameter(atom::electronegativity)(molecule.atoms()[i]);
        chi0[i] = chi;
        q[i] = - chi;
    }

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n_int, n_int, m_int, 1.0, B, n_int, tmp, n_int, 1.0, res, n_int);

    int info = LAPACKE_dsysv(LAPACK_ROW_MAJOR, 'U', n_int, 1, res, n_int, ipiv, chi0, 1);
    if (info) {
        throw std::runtime_error("Cannot solve linear system");
    }

    for(size_t i = 0; i < n; i++) {
        q[i] += chi0[i];
    }

    mkl_free(W);
    mkl_free(B);
    mkl_free(tmp);
    mkl_free(res);
    mkl_free(chi0);
    mkl_free(ipiv);

    return q;
}
