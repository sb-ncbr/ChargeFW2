//
// Created by krab1k on 23/10/18.
//

#include <utility>
#include "atom.h"
#include "../element.h"
#include "../periodic_table.h"


Atom::Atom(int idx, const Element *element, double x, double y, double z) {
    index_ = idx;
    element_ = element;
    pos_[0] = x;
    pos_[1] = y;
    pos_[2] = z;
}
