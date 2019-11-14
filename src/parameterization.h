#include <utility>

//
// Created by krab1k on 20.11.18.
//

#pragma once

#include <string>

#include "structures/molecule_set.h"
#include "method.h"
#include "charges.h"
#include "parameters.h"


class Parameterization {
    const MoleculeSet &set_;
    Method *method_;
    const Charges &reference_charges_;
    std::unique_ptr<Parameters> parameters_;
    std::string charge_output_dir_;
    std::string parameters_output_file_;

public:
    Parameterization(const MoleculeSet &ms, Method *method, const Charges &reference_charges,
                     const std::string &charge_output_file, const std::string &parameters_output_file);

    void parametrize();

    [[nodiscard]] const MoleculeSet &set() const { return set_; }

    [[nodiscard]] Method *method() const { return method_; }

    [[nodiscard]] const Charges &reference_charges() const { return reference_charges_; }
};
