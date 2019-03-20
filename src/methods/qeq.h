//
// Created by krab1k on 02/01/19.
//

#pragma once

#include <string>
#include <vector>
#include <boost/config.hpp>

#include "../structures/atom.h"
#include "../structures/molecule.h"
#include "../method.h"


class QEq : public EEMethod {
    enum atom{electronegativity, hardness};
    double overlap_term(const Atom &atom_i, const Atom &atom_j, const std::string &type) const;
    std::vector<double> solve_system(const std::vector<const Atom *> &atoms, double total_charge) const override;
public:
    explicit QEq() : EEMethod("QEq", {}, {"electronegativity", "hardness"}, {},
            {
                    {"overlap_term", {"overlap_term", "Overlap term", "str", "Louwen-Vogt",
                                             {"Nishimoto-Mataga",
                                              "Nishimoto-Mataga-Weiss",
                                              "Ohno",
                                              "Ohno-Klopman",
                                              "DasGupta-Huzinaga",
                                              "Louwen-Vogt"}}}
            }) {}

    virtual ~QEq() = default;
};

extern "C" BOOST_SYMBOL_EXPORT QEq method;
QEq method;
