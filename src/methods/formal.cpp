#include "formal.h"
#include "../method_registry.h"


[[maybe_unused]] const bool Formal_registered_ =
    (MethodRegistry::register_factory("formal", &make_method<Formal>), true);

std::vector<double> Formal::calculate_charges(const Molecule &molecule) const {
    std::vector<double> res;
    res.reserve(molecule.atoms().size());
    for (const auto &atom: molecule.atoms()) {
        res.push_back(atom.formal_charge());
    }
    return res;
}
