//
// Created by krab1k on 24.1.19.
//

#pragma once

#include "reader.h"
#include "writer.h"
#include "../charges.h"


class Mol2 : public Reader, public Writer {

public:
    MoleculeSet read_file(const std::string &filename, bool) override;

    void save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) override;
};
