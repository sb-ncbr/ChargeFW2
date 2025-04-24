#pragma once

#include <Eigen/Core>
#include <string>
#include <vector>

#include "../structures/atom.h"
#include "../structures/molecule.h"
#include "../method.h"


class QEq : public EEMethod {
    inline static const MethodMetadata METADATA = {
        .name = "QEq",
        .internal_name = "qeq",
        .full_name = "Charge Equilibration",
        .publication = "10.1021/j100161a070",
        .type = "3D",
        .priority = 170
    };

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

    [[nodiscard]] const MethodMetadata& metadata() const override {
        return METADATA;
    }

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
