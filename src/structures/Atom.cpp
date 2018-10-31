//
// Created by krab1k on 23/10/18.
//

#include <utility>
#include "Atom.h"
#include "../Element.h"

#include <iostream>

using std::cout;
using std::endl;


Atom::Atom(int idx, const Element *element, double x, double y, double z) {
    index_ = idx;
    element_ = element;
    pos_[0] = x;
    pos_[1] = y;
    pos_[2] = z;
}

bool Atom::operator==(const Atom &other) const {
    if (index_ != other.index_)
        return false;

    if (element_ != other.element_)
        return false;

    return !(pos_ != other.pos_);

}

std::ostream &operator<<(std::ostream &str, const Atom &atom) {
    return str << "Atom " << atom.element_->symbol().toStdString() << " Idx: " << atom.index_;
}
