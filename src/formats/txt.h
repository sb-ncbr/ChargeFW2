#pragma once

#include <string>

#include "writer.h"
#include "../charges.h"


class TXT final: public Writer {
public:
    void save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) override;
};
