#pragma once

#include <string>

#include "../charges.h"
#include "../structures/molecule_set.h"

void save_charges(const MoleculeSet& ms, const Charges& charges, const std::string& filename);
