//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <utility>

#include "Atom.h"

class Bond {
    Atom first_;
    Atom second_;
    int order_{};
public:
    Bond() = default;
    Bond(Atom &atom1, Atom &atom2, int order);
};

