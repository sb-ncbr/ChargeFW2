//
// Created by krab1k on 02/01/19.
//

#pragma once

#include <Eigen/Core>
#include <string>
#include <vector>
#include <boost/config.hpp>

#include "../structures/atom.h"
#include "../structures/molecule.h"
#include "../method.h"


class QEq : public EEMethod {
    enum atom{electronegativity, hardness};
    [[nodiscard]] double overlap_term(const Atom &atom_i, const Atom &atom_j, const std::string &type) const;

    [[nodiscard]] Eigen::VectorXd EE_system(const std::vector<const Atom *> &atoms, double total_charge) const;

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

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
