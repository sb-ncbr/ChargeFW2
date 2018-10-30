//
// Created by krab1k on 29/10/18.
//

#include <iostream>

#include "Bond.h"

using std::cout;
using std::endl;

Bond::Bond(Atom &atom1, Atom &atom2, int order) {
    first_ = atom1;
    second_ = atom2;
    order_ = order;
}

void Bond::print() const {
    cout << *this << endl;
}

std::ostream &operator<<(std::ostream &str, const Bond &bond) {
    return str << "Bond: (" << bond.first_ << ", " << bond.second_ << ")";
}
