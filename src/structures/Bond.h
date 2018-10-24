//
// Created by krab1k on 24/10/18.
//

#pragma once

#include "Atom.h"

class Bond {
    const Atom &first_;
    const Atom &second_;
    int order_;
public:
    Bond(const Atom &a1, const Atom &a2, int order) : first_{a1}, second_{a2}, order_{order} {};
};


