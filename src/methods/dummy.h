//
// Created by krab1k on 8.11.18.
//

#pragma once

#include "../method.h"


class Dummy : public Method {
public:
    explicit Dummy() : Method("Dummy", {}, {}, {}, {}) {}

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
