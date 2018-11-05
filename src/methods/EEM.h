//
// Created by krab1k on 31/10/18.
//

#pragma once

#include "../structures/Molecule.h"
#include "../Method.h"

class EEM : public Method {
public:
    void calculate_charges(const Molecule &molecule) override;
};
