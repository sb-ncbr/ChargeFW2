//
// Created by krab1k on 8.11.18.
//

#include "Dummy.h"

std::vector<double> Dummy::calculate_charges(const Molecule &molecule) {
    return std::vector<double>(molecule.atoms().size(), 0);
}
