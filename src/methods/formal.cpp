#include "formal.h"
#include <optional>

CHARGEFW2_METHOD(Formal)

namespace {
    const MethodMetadata FORMAL_METADATA = {
        .internal_name = "formal",
        .full_name = "Formal charges",
        .publication = std::nullopt,
        .type = "other",
        .priority = 10
    };
};

const MethodMetadata& Formal::get_metadata() const {
    return FORMAL_METADATA;
};

std::vector<double> Formal::calculate_charges(const Molecule &molecule) const {
    std::vector<double> res;
    res.reserve(molecule.atoms().size());
    for (const auto &atom: molecule.atoms()) {
        res.push_back(atom.formal_charge());
    }
    return res;
}
