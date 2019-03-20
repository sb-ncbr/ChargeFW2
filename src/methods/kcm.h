//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"


class KCM : public Method {
    enum atom{electronegativity, hardness};
public:
    explicit KCM() : Method("KCM", {}, {"electronegativity", "hardness"}, {}, {}) {}

    virtual ~KCM() = default;

    std::vector<double> calculate_charges(const Molecule &molecule) const override;

    bool is_suitable_for_large_molecule() const override { return false; }
};

extern "C" BOOST_SYMBOL_EXPORT KCM method;
KCM method;
