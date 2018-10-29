//
// Created by krab1k on 29/10/18.
//

#include "Bond.h"

Bond::Bond(Atom &atom1, Atom &atom2, int order) {
    first_ = atom1;
    second_ = atom2;
    order_ = order;
}