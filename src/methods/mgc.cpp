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

    double product = 1;
    for (size_t i = 0; i < n; i++) {
        #define IDX(i, j) ((i) * n + (j))
        auto &atom_i = molecule.atoms()[i];
        S[IDX(i, i)] = 1 + molecule.degree(atom_i);
        for(size_t j = i + 1; j < n; j++) {
            auto &atom_j = molecule.atoms()[j];
            S[IDX(i, j)] = -molecule.bond_order(atom_i, atom_j);
        }
        X0[i] = atom_i.element().electronegativity();
        product *= X0[i];
    }

    auto n_int = static_cast<int>(n);
    int info = LAPACKE_dsysv(LAPACK_ROW_MAJOR, 'U', n_int, 1, S, n_int, ipiv, X0, 1);
    if(info) {
        throw std::runtime_error("Cannot solve linear system");
    }

    std::vector<double> results;
    results.assign(X0, X0 + n);

    product = pow(product, 1.0 / n);

    for (size_t i = 0; i < n; i++) {
        results[i] -= molecule.atoms()[i].element().electronegativity();
        results[i] /= product;
    }

    mkl_free(S);
    mkl_free(X0);
    mkl_free(ipiv);

    return results;
}
