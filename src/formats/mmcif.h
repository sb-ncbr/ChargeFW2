//
// Created by krab1k on 28.1.19.
//


#pragma once

#include "reader.h"


class mmCIF: public Reader {

public:
    MoleculeSet read_file(const std::string &filename, bool read_hetatms, bool ignore_water) override;
};
