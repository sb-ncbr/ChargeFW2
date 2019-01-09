//
// Created by krab1k on 02/01/19.
//

#include <string>
#include <cmath>
#include <mkl_lapacke.h>
#include <mkl.h>

#include "../structures/atom.h"
#include "qeq.h"
#include "config.h"
#include "../geometry.h"
#include "../parameters.h"

#define IDX(i, j) ((i) * m + (j))


double QEq::overlap_term(const Atom &atom_i, const Atom &atom_j, std::string type) {
    auto Ji = parameters_->atom()->parameter(atom::hardness)(atom_i);
    auto Jj = parameters_->atom()->parameter(atom::hardness)(atom_j);
    auto Rij = distance(atom_i, atom_j);
    if (type == "Nishimoto-Mataga") {
        return 1 / (Rij + 2 / (Ji + Jj));
    } else if (type == "Nishimoto-Mataga-Weiss") {
        const double f = 1.2;
        return f / (Rij + (2 * f) / (Ji + Jj));
    } else if (type == "Ohno") {
        return 1 / std::sqrt(Rij * Rij + std::pow(2 / (Ji + Jj), 2));
    } else if (type == "Ohno-Klopman") {
        return 1 / std::sqrt(Rij * Rij + std::pow(1 / (2 * Ji) + 1 / (2 * Jj), 2));
    } else if (type == "DasGupta-Huzinaga") {
        const double k = 0.4;
        return 1 / (Rij + 1 / (Ji / 2 * exp(k * Rij) + Jj / 2 * exp(k * Rij)));
    } else if (type == "Louwen-Vogt") {
        const double gamma = (Ji + Jj) / 2;
        return 1 / std::cbrt(1 / std::pow(gamma, 3) + std::pow(Rij, 3));
    } else {
        std::cerr << "Unknown overlap term type: " << type << std::endl;
        exit(EXIT_PARAMETER_ERROR);
    }
}


std::vector<double> QEq::calculate_charges(const Molecule &molecule) {

    size_t n = molecule.atoms().size();
    size_t m = n + 1;

    auto *A = (double *) mkl_malloc(m * m * sizeof(double), 64);
    auto *b = (double *) mkl_malloc(m * sizeof(double), 64);
    auto *ipiv = (MKL_INT *) mkl_malloc(m * sizeof(MKL_INT), 64);

    for (size_t i = 0; i < n; i++) {
        const auto &atom_i = molecule.atoms()[i];
        A[IDX(i, i)] = parameters_->atom()->parameter(atom::hardness)(atom_i);
        b[i] = - parameters_->atom()->parameter(atom::electronegativity)(atom_i);
        for (size_t j = i + 1; j < n; j++) {
            const auto &atom_j = molecule.atoms()[j];
            auto type = get_option_value<std::string>("overlap_term");
            A[IDX(i, j)] = overlap_term(atom_i, atom_j, type);
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