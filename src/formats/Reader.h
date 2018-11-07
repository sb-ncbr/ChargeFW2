//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <string>
#include "../structures/MoleculeSet.h"

class Reader {
public:
    virtual MoleculeSet read_file(const std::string &filename) = 0;
};


