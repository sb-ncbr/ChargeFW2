#include "dummy.h"

CHARGEFW2_METHOD(Dummy)


std::vector<double> Dummy::calculate_charges(const Molecule &molecule) const {
    return std::vector<double>(molecule.atoms().size(), 0);
}
