//
// Created by krab1k on 8.11.18.
//

#include "veem.h"

std::vector<double> VEEM::calculate_charges(const Molecule &molecule) {
    size_t n = molecule.atoms().size();
    std::vector<double> q(n, 0);

    double num = 0;
    double den = 0;

    for(size_t i = 0; i < n; i++) {
        auto &atom = molecule.atoms()[i];
        num += atom.element().electronegativity() * atom.element().valence_electron_count();
        den += atom.element().valence_electron_count();
    }

    double eq_en = num / den;

    for(size_t i = 0; i < n; i++) {
        auto &atom = molecule.atoms()[i];
        q[i] = atom.element().valence_electron_count() * (eq_en - atom.element().electronegativity()) / eq_en;
    }

    return q;
}
