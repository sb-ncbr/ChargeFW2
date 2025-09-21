#include "dummy.h"

#include "../method_registry.h"


[[maybe_unused]] const bool Dummy_registered_ =
    (MethodRegistry::register_factory("dummy", &make_method<Dummy>), true);
std::vector<double> Dummy::calculate_charges(const Molecule &molecule) const {
    return std::vector<double>(molecule.atoms().size(), 0);
}
