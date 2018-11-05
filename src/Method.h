//
// Created by krab1k on 31/10/18.
//

#pragma once

#include "structures/Molecule.h"

class Method {

public:
    virtual void calculate_charges(const Molecule &molecule) = 0;
};


