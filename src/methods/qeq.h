#pragma once

#include <Eigen/Core>
#include <string>
#include <vector>

#include "../structures/atom.h"
#include "../structures/molecule.h"
#include "../method.h"


class QEq : public EEMethod {
    enum atom{electronegativity, hardness};
    [[nodiscard]] double overlap_term(const Atom &atom_i, const Atom &atom_j, const std::string &type) const;

    [[nodiscard]] Eigen::VectorXd EE_system(const std::vector<const Atom *> &atoms, double total_charge) const;

public:
    explicit QEq() : EEMethod({}, {"electronegativity", "hardness"}, {},
            {
                    {"overlap_term", {"overlap_term", "Overlap term", "str", "Louwen-Vogt",
                                             {"Nishimoto-Mataga",
                                              "Nishimoto-Mataga-Weiss",
                                              "Ohno",
                                              "Ohno-Klopman",
                                              "DasGupta-Huzinaga",
                                              "Louwen-Vogt"}}}
            }) {}

    [[nodiscard]] const MethodMetadata& metadata() const override;

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
