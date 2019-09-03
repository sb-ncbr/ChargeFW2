//
// Created by krab1k on 12/11/18.
//

#pragma once

#include <boost/config.hpp>

#include "../method.h"


class Formal : public Method {
public:
    explicit Formal() : Method("Formal", {}, {}, {}, {}) {}

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

    virtual ~Formal() = default;
};

extern "C" BOOST_SYMBOL_EXPORT Formal method;
Formal method;
