#include "veem.h"

CHARGEFW2_METHOD(VEEM)

namespace {
    const MethodMetadata VEEM_METADATA = {
        .internal_name = "veem",
        .full_name = "Valence Electrons Equalization Method",
        .publication = "10.1088/1674-0068/24/01/31-39",
        .type = "2D",
        .priority = 20
    };
};

const MethodMetadata& VEEM::get_metadata() const {
    return VEEM_METADATA;
};

std::vector<double> VEEM::calculate_charges(const Molecule &molecule) const {
    size_t n = molecule.atoms().size();
    std::vector<double> q(n, 0);

    double num = 0;
    double den = 0;

    for (size_t i = 0; i < n; i++) {
        auto &atom = molecule.atoms()[i];
        num += atom.element().electronegativity() * atom.element().valence_electron_count();
        den += atom.element().valence_electron_count();
    }

    double eq_en = num / den;

    for (size_t i = 0; i < n; i++) {
        auto &atom = molecule.atoms()[i];
        q[i] = atom.element().valence_electron_count() * (eq_en - atom.element().electronegativity()) / eq_en;
    }

    return q;
}


bool VEEM::is_suitable_for_molecule(const Molecule &molecule) const {
    try {
        for (const auto &atom: molecule.atoms()) {
            (void) atom.element().valence_electron_count();
        }
        return true;
    }
    catch (const std::exception &) {
        return false;
    }
}
