#include "dummy.h"
#include <optional>

CHARGEFW2_METHOD(Dummy)

namespace {
    const MethodMetadata DUMMY_METADATA = {
        .name = "Dummy",
        .internal_name = "dummy",
        .full_name = "Dummy Method",
        .publication = std::nullopt,
        .type = "other",
        .priority = 0
    };
};

const MethodMetadata& Dummy::metadata() const {
    return DUMMY_METADATA;
};

std::vector<double> Dummy::calculate_charges(const Molecule &molecule) const {
    return std::vector<double>(molecule.atoms().size(), 0);
}
