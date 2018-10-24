//
// Created by krab1k on 23/10/18.
//

#include <string>

#include "Atom.h"

Atom::Atom(int idx, std::string symbol, double x, double y, double z) {
    index_ = idx;
    symbol_ = symbol;
    pos_[0] = x;
    pos_[1] = y;
    pos_[2] = z;
}
