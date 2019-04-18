//
// Created by krab1k on 31/10/18.
//

#include <vector>
#include <cmath>
#include <mkl_lapacke.h>
#include <mkl.h>

#include "mgc.h"
#include "../structures/molecule.h"


std::vector<double> MGC::calculate_charges(const Molecule &molecule) const {

    size_t n = molecule.atoms().size();

    auto *S = static_cast<double *>(mkl_malloc(n * n * sizeof(double), 64));
    auto *X0 = static_cast<double *>(mkl_malloc(n * sizeof(double), 64));
    auto *ipiv = static_cast<int *>(mkl_malloc(n * sizeof(int), 64));

    double log_sum = 0;

    #define IDX(i, j) ((i) * n + (j))
    for (const auto &atom: molecule.atoms()) {
        auto i = atom.index();
        S[IDX(i, i)] = 1;
        X0[i] = atom.element().electronegativity();
        log_sum += log(X0[i]);
    }

    for (const auto &bond: molecule.bonds()) {
        auto i1 = bond.first().index();
        auto i2 = bond.second().index();
        auto order = bond.order();
        S[IDX(i1, i1)] += order;
        S[IDX(i2, i2)] += order;
        S[IDX(i1, i2)] -= order;
        S[IDX(i2, i1)] -= order;
    }

    auto n_int = static_cast<int>(n);
    int info = LAPACKE_dsysv(LAPACK_ROW_MAJOR, 'U', n_int, 1, S, n_int, ipiv, X0, 1);
    if(info) {
        throw std::runtime_error("Cannot solve linear system");
    }

    std::vector<double> results;
    results.assign(X0, X0 + n);

    double avg_chi = exp(log_sum / n);

    for (size_t i = 0; i < n; i++) {
        results[i] -= molecule.atoms()[i].element().electronegativity();
        results[i] /= avg_chi;
    }

    mkl_free(S);
    mkl_free(X0);
    mkl_free(ipiv);

    return results;
}
