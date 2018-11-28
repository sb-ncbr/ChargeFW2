//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <string>
#include "../structures/molecule_set.h"

class reader {
public:
    virtual MoleculeSet read_file(const std::string &filename) = 0;
};


