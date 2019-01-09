//
// Created by krab1k on 29/10/18.
//

#include "../element.h"
#include "bond.h"


std::array<double, 3> Bond::get_center(bool weighted) const {
    auto &p1 = first_->pos();
    auto &p2 = second_->pos();

    std::array<double, 3> pos{};

    if (weighted) {
        double c1 = first_->element().covalent_radius();
        double c2 = second_->element().covalent_radius();
        pos = {(c1 * p1[0] + c2 * p2[0]) / (c1 + c2), (c1 * p1[1] + c2 * p2[1]) / (c1 + c2),
               (c1 * p1[2] + c2 * p2[2]) / (c1 + c2)};

    } else {
        pos = {(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2};
    }

    return pos;
}
