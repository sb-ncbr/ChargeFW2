//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"

class EEM : public Method {
    enum common {kappa};
    enum atom {A, B};

    std::vector<double> solve_EEM_system(const std::vector<const Atom *> &atoms, double total_charge) const;
public:
    explicit EEM() : Method("EEM", {"kappa"}, {"A", "B"}, {},
            {
                    {"type", {"type", "Type of a solver", "str", "full", {"full", "cutoff"}}},
                    {"radius", {"radius", "Radius for cutoff", "double", "8", {}}}
            }) {}

    std::vector<double> calculate_charges(const Molecule &molecule) const override;

    std::vector<RequiredFeatures> get_requirements() const override {
        if (get_option_value<std::string>("type") == "cutoff") {
            return {RequiredFeatures::DISTANCE_TREE};
        } else {
            return std::vector<RequiredFeatures>();
        }
    }
};

extern "C" BOOST_SYMBOL_EXPORT EEM method;
EEM method;