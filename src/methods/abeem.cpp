//
// Created by krab1k on 31/10/18.
//

#include <vector>
#include <cmath>
#include <mkl_lapacke.h>
#include <mkl.h>

#include "abeem.h"
#include "../parameters.h"
#include "../geometry.h"


#define IDX(i, j) ((i) * mn + (j))

std::vector<double> ABEEM::calculate_charges(const Molecule &molecule) {

    size_t n = molecule.atoms().size();
    size_t m = molecule.bonds().size();
    size_t mn = n + m + 1;

    auto *A = (double *) mkl_malloc(mn * mn * sizeof(double), 64);
    auto *b = (double *) mkl_malloc(mn * sizeof(double), 64);
    auto *ipiv = (MKL_INT *) mkl_malloc(mn * sizeof(MKL_INT), 64);

    const double k = parameters_->common()->parameter(common::k);

    // atom-atom part
    for (size_t i = 0; i < n; i++) {
        const auto &atom_i = molecule.atoms()[i];
        A[IDX(i, i)] = parameters_->atom()->parameter(atom::b)(atom_i);
        b[i] = -parameters_->atom()->parameter(atom::a)(atom_i);
        for (size_t j = i + 1; j < n; j++) {
            const auto &atom_j = molecule.atoms()[j];
            double off = k / distance(atom_i, atom_j);
            A[IDX(i, j)] = off;
            A[IDX(j, i)] = off;
        }
    }

    // atom-bond part
    for (size_t i = 0; i < n; i++) {
        const auto &atom = molecule.atoms()[i];
        for (size_t j = 0; j < m; j++) {
            const auto &bond = molecule.bonds()[j];
            if (bond.hasAtom(atom)) {
                A[IDX(i, n + j)] = parameters_->atom()->parameter(atom::c)(atom);
            } else {
                A[IDX(i, n + j)] = k / distance(atom, bond, true);
            }
        }

    }

    // bond-atom part
    for (size_t i = 0; i < m; i++) {
        const auto &bond = molecule.bonds()[i];
        b[n + i] = -parameters_->bond()->parameter(bond::A)(bond);
        for (size_t j = 0; j < n; j++) {
            const auto &atom = molecule.atoms()[j];
            if (bond.hasAtom(atom)) {
                if (bond.first() == atom) {
                    A[IDX(n + i, j)] = parameters_->bond()->parameter(bond::D)(bond);
                } else {
                    A[IDX(n + i, j)] = parameters_->bond()->parameter(bond::C)(bond);
                }
            } else {
                A[IDX(n + i, j)] = k / distance(atom, bond, true);
            }
        }
    }

    // bond-bond part
    for (size_t i = 0; i < m; i++) {
        const auto &bond_i = molecule.bonds()[i];
        A[IDX(n + i, n + i)] = parameters_->bond()->parameter(bond::B)(bond_i);
        for (size_t j = i + 1; j < m; j++) {
            const auto &bond_j = molecule.bonds()[j];
            double off = k / distance(bond_i, bond_j, true);
            A[IDX(n + i, n + j)] = off;
            A[IDX(n + j, n + i)] = off;
        }
    }

    for (size_t i = 0; i < n + m; i++) {
        A[IDX(i, n + m)] = 1;
        A[IDX(n + m, i)] = 1;
    }

    A[IDX(n + m, n + m)] = 0;
    b[n + m] = molecule.total_charge();

    MKL_INT info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, mn, 1, A, mn, ipiv, b, 1);
    if (info) {
        throw std::runtime_error("Cannot solve linear system");
    }

    // Redistribute the bond charges to the corresponding atoms
    for(size_t i = 0; i < m; i++) {
        const auto &bond = molecule.bonds()[i];
        b[bond.first().index()] += 0.5 * b[n + i];
        b[bond.second().index()] += 0.5 * b[n + i];
    }

    std::vector<double> results;
    results.assign(b, b + n);

    mkl_free(A);
    mkl_free(b);
    mkl_free(ipiv);

    return results;
}