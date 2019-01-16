//
// Created by krab1k on 31/10/18.
//

#include <vector>
#include <cmath>
#include <string>
#include <mkl_lapacke.h>
#include <mkl.h>

#include "eem.h"
#include "../parameters.h"
#include "../geometry.h"

#define IDX(i, j) ((i) * m + (j))


std::vector<double> EEM::solve_EEM_system(const std::vector<const Atom *> atoms, double total_charge) const {

    size_t n = atoms.size();
    size_t m = n + 1;

    auto *A = (double *) mkl_malloc(m * m * sizeof(double), 64);
    auto *b = (double *) mkl_malloc(m * sizeof(double), 64);
    auto *ipiv = (MKL_INT *) mkl_malloc(m * sizeof(MKL_INT), 64);

    for (size_t i = 0; i < n; i++) {
        const auto &atom_i = *atoms[i];
        A[IDX(i, i)] = parameters_->atom()->parameter(atom::B)(atom_i);
        b[i] = -parameters_->atom()->parameter(atom::A)(atom_i);
        for (size_t j = i + 1; j < n; j++) {
            const auto &atom_j = *atoms[j];
            A[IDX(i, j)] = parameters_->common()->parameter(common::kappa) / distance(atom_i, atom_j);
        }
    }

    for (size_t i = 0; i < n; i++) {
        A[IDX(i, n)] = 1;
    }

    A[IDX(n, n)] = 0;
    b[n] = total_charge;

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


std::vector<double> EEM::calculate_charges(const Molecule &molecule) {

    auto method = get_option_value<std::string>("type");
    std::vector<const Atom *> fragment_atoms;
    if (method == "full") {
        for (const auto &atom: molecule.atoms()) {
            fragment_atoms.push_back(&atom);
        }

        return solve_EEM_system(fragment_atoms, molecule.total_charge());

    } else /* method == "cutoff" */ {
        std::vector<double> results;
        for (const auto &atom: molecule.atoms()) {
            fragment_atoms = molecule.get_close_atoms(atom, get_option_value<double>("radius"));
            auto res = solve_EEM_system(fragment_atoms,
                                        static_cast<double>(molecule.total_charge()) * fragment_atoms.size() /
                                        molecule.atoms().size());
            results.push_back(res[0]);
        }

        double correction = molecule.total_charge();
        for (auto val: results) {
            correction -= val;
        }

        correction /= molecule.atoms().size();

        for (auto &val: results) {
            val += correction;
        }

        return results;
    }
}
