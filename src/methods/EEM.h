//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>

#include "../structures/Molecule.h"
#include "../Method.h"

class EEM : public Method {
    enum common{kappa};
    enum atom{A, B};
public:
    explicit EEM(const Parameters *parameters) : Method({"kappa"}, {"A", "B"}, {}, parameters) {}

    std::vector<double> calculate_charges(const Molecule &molecule) override;
};
