//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <string>

#include "reader.h"

class SDF: public reader {

public:
    MoleculeSet read_file(const std::string &filename) override;
};