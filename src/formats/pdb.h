#pragma once

#include "reader.h"

class PDB final : public Reader {

public:
    PDB();
    MoleculeSet read_file(const std::string &filename) override;
};
