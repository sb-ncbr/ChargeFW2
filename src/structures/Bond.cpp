//
// Created by krab1k on 29/10/18.
//

#include <iostream>

#include "Bond.h"


std::ostream &operator<<(std::ostream &str, const Bond &bond) {
    return str << "Bond: (" << *bond.first_ << ", " << *bond.second_ << "): " << bond.order_;
}
