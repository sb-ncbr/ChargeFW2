//
// Created by krab1k on 02/01/19.
//

#include <cmath>

#include "structures/atom.h"
#include "geometry.h"

double distance(const Atom &atom1, const Atom &atom2) {
    double dx = atom1.pos()[0] - atom2.pos()[0];
    double dy = atom1.pos()[1] - atom2.pos()[1];
    double dz = atom1.pos()[2] - atom2.pos()[2];

    return std::sqrt(dx * dx + dy * dy + dz * dz);
}