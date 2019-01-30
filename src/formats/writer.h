//
// Created by krab1k on 30.1.19.
//

#pragma once

#include <string>

#include "../charges.h"
#include "../structures/molecule_set.h"


class Writer {

public:
    virtual void save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) = 0;

    virtual ~Writer();
};
