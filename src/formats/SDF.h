//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <QString>

#include "Reader.h"

class SDF: public Reader {

public:
    MoleculeSet read_file(const QString &filename) override;
};


