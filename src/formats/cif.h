//
// Created by danny305 on May 3rd 2021.
//

#pragma once

#include <string>
#include <gemmi/cif.hpp>

#include "writer.h"
#include "../structures/molecule_set.h"
#include "../charges.h"


class CIF : public Writer {
public:
    void save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) override;
};
