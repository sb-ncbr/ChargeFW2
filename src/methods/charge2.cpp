//
// Created by krab1k on 13.11.18.
//

#include <cmath>

#include "charge2.h"
#include "../structures/molecule.h"
#include "../structures/bond.h"

const int n_iters = 10;


std::vector<double> Charge2::calculate_charges(const Molecule &molecule) {

    const size_t n = molecule.atoms().size();
    std::vector<double> q(n, 0);
    for (int i = 0; i < n_iters; i++) {
        for (auto &atom: molecule.atoms()) {
            double alpha_charge = 0.0;
            for (auto &bonded: molecule.k_bond_distance(atom, 1)) {
                double a;
                if (atom.element().period() == 2 and bonded->element().period() == 2) {
                    a = parameters_->common()->parameter(common::a1);
                } else if (atom.element().symbol() == "H" or bonded->element().symbol() == "H") {
                    a = parameters_->common()->parameter(common::a2);
                } else {
                    a = parameters_->common()->parameter(common::a3);
                }
                alpha_charge += (parameters_->atom()->parameter(atom::chi)(*bonded) -
                                 parameters_->atom()->parameter(atom::chi)(*bonded)) / a;
            }
            double beta_charge = 0.0;
            double P = parameters_->atom()->parameter(atom::P0)(atom) * (1 + parameters_->common()->parameter(
                    common::alpha) * (parameters_->atom()->parameter(atom::q0)(atom) - q[atom.index()]));
            for (auto &bonded: molecule.k_bond_distance(atom, 2)) {
                // 7.17 should be replaced with something like:
                // parameters->atom()->parameter(atom::chi)(hydrogen)
                beta_charge += (parameters_->atom()->parameter(atom::chi)(*bonded) - 7.17) * P /
                               parameters_->common()->parameter(common::b);
            }
            double gamma_charge = 0.0;
            for (auto &bonded: molecule.k_bond_distance(atom, 3)) {
                beta_charge += (parameters_->atom()->parameter(atom::chi)(*bonded) - 7.17) * P /
                               parameters_->common()->parameter(common::b) /
                               parameters_->common()->parameter(common::c);
            }

            q[atom.index()] = alpha_charge + beta_charge + gamma_charge;
        }

    }

    return q;
}