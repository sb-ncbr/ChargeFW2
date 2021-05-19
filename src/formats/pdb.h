//
// Created by krab1k on 24.1.19.
//

#pragma once

#include "reader.h"

class PDB : public Reader {

public:
    PDB();
    MoleculeSet read_file(const std::string &filename) override;
};
