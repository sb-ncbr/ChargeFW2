//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <QString>
#include "../structures/MoleculeSet.h"

class Reader {
public:
    virtual MoleculeSet read_file(const QString &filename) = 0;
};


