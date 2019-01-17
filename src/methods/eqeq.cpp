//
// Created by krab1k on 31/10/18.
//

#include <vector>
#include <cmath>
#include <mkl_lapacke.h>
#include <mkl.h>

#include "eqeq.h"
#include "../parameters.h"
#include "../geometry.h"


#define IDX(i, j) ((i) * m + (j))

std::vector<double> EQeq::calculate_charges(const Molecule &molecule) const {

    size_t n = molecule.atoms().size();
    size_t m = n + 1;

    const double lambda = 1.2;
    const double k = 14.4;
    double H_electron_affinity = -2.0; // Exception for hydrogen mentioned in the article

    auto *A = (double *) mkl_malloc(m * m * sizeof(double), 64);
    auto *b = (double *) mkl_malloc(m * sizeof(double), 64);
    auto *ipiv = (MKL_INT *) mkl_malloc(m * sizeof(MKL_INT), 64);

    std::vector<double> X(n, 0);
    std::vector<double> J(n, 0);

    for (size_t i = 0; i < n; i++) {
        const auto &atom_i = molecule.atoms()[i];
        if (atom_i.element().symbol() == "H") {
            X[i] = (atom_i.element().ionization_potential() + H_electron_affinity) / 2;
            J[i] = atom_i.element().ionization_potential() - H_electron_affinity;
        } else {
            X[i] = (atom_i.element().ionization_potential() + atom_i.element().electron_affinity()) / 2;
            J[i] = atom_i.element().ionization_potential() - atom_i.element().electron_affinity();
        }
    }

    for (size_t i = 0; i < n; i++) {
        const auto &atom_i = molecule.atoms()[i];
        A[IDX(i, i)] = J[i];
        b[i] = -X[i];
        for (size_t j = i + 1; j < n; j++) {
            const auto &atom_j = molecule.atoms()[j];
            double a = std::sqrt(J[i] * J[j]) / k;
            double Rij = distance(atom_i, atom_j);
            double overlap = std::exp(-a * a * Rij * Rij) * (2 * a - a * a * Rij - 1 / Rij);
            A[IDX(i, j)] = lambda * k / 2 * (1 / Rij + overlap);
        }
    }

    for (size_t i = 0; i < n; i++) {
        A[IDX(i, n)] = 1;
    }

    A[IDX(n, n)] = 0;
    b[n] = molecule.total_charge();

    MKL_INT info = LAPACKE_dsysv(LAPACK_ROW_MAJOR, 'U', m, 1, A, m, ipiv, b, 1);
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