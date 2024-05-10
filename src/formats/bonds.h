#pragma once

#include <memory>
#include <vector>

#include "../structures/atom.h"
#include "../structures/bond.h"


std::unique_ptr<std::vector<Bond>> get_bonds(std::unique_ptr<std::vector<Atom>> &atoms);
