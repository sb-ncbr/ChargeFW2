//
// Created by krab1k on 02/01/19.
//

#include <cmath>

#include "structures/atom.h"
#include "structures/bond.h"
#include "geometry.h"


double distance(const Atom &atom1, const Atom &atom2) {
    double dx = atom1.pos()[0] - atom2.pos()[0];
    double dy = atom1.pos()[1] - atom2.pos()[1];
    double dz = atom1.pos()[2] - atom2.pos()[2];

    return std::sqrt(dx * dx + dy * dy + dz * dz);
}


double distance(const Atom &atom, const Bond &bond, bool weighted) {
    std::array<double, 3> bc = bond.get_center(weighted);
    double dx = atom.pos()[0] - bc[0];
    double dy = atom.pos()[1] - bc[1];
    double dz = atom.pos()[2] - bc[2];

    return std::sqrt(dx * dx + dy * dy + dz * dz);
}


double distance(const Bond &bond1, const Bond &bond2, bool weighted) {
    std::array<double, 3> bc1 = bond1.get_center(weighted);
    std::array<double, 3> bc2 = bond2.get_center(weighted);
    double dx = bc1[0] - bc2[0];
    double dy = bc1[1] - bc2[1];
    double dz = bc1[2] - bc2[2];

    return std::sqrt(dx * dx + dy * dy + dz * dz);
}
