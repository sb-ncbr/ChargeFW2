//
// Created by krab1k on 13.11.18.
//

#include <cmath>

#include "gdac.h"
#include "../structures/molecule.h"
#include "../structures/bond.h"
#include "../geometry.h"
#include "../parameters.h"

CHARGEFW2_METHOD(GDAC)


std::vector<double> GDAC::calculate_charges(const Molecule &molecule) const {

    const size_t n = molecule.atoms().size();
    std::vector<double> q(n, 0);
    std::vector<double> chi(n, 0);

    for (int alpha = 1; alpha < get_option_value<int>("iters"); alpha++) {
        for (size_t i = 0; i < n; i++) {
            const auto &atom = molecule.atoms()[i];
            chi[i] = parameters_->atom()->parameter(atom::B)(atom) * q[i] +
                     parameters_->atom()->parameter(atom::A)(atom);
        }

        for (const auto &bond: molecule.bonds()) {
            const Atom *atom1 = &bond.first();
            const Atom *atom2 = &bond.second();

            double chi1 = chi[atom1->index()];
            double chi2 = chi[atom2->index()];

            if (chi1 > chi2) {
                atom1 = &bond.second();
                atom2 = &bond.first();
                std::swap(chi1, chi2);
            }

            double d = parameters_->atom()->parameter(atom::A)(*atom1) +
                       parameters_->atom()->parameter(atom::B)(*atom1);

            double f = 1 - (distance(*atom1, *atom2) / (atom1->element().vdw_radius() + atom2->element().vdw_radius()));
            double diff = pow(f, alpha) * (chi2 - chi1) / d;
            q[atom1->index()] += diff;
            q[atom2->index()] -= diff;
        }
    }

    return q;
}
