#include <utility>

//
// Created by krab1k on 20.11.18.
//

#pragma once

#include <boost/shared_ptr.hpp>
#include <string>

#include "structures/molecule_set.h"
#include "method.h"
#include "charges.h"
#include "parameters.h"


class Parameterization {
    const MoleculeSet &set_;
    boost::shared_ptr<Method> method_;
    const Charges &reference_charges_;
    std::unique_ptr<Parameters> parameters_;
    std::string charge_output_file_;


public:
    Parameterization(const MoleculeSet &ms, boost::shared_ptr<Method> method, const Charges &reference_charges,
                     const std::string &charge_output_file);

    void parametrize();

    const MoleculeSet &set() const { return set_; }

    boost::shared_ptr<Method> method() { return method_; }

    const Charges &reference_charges() const { return reference_charges_; }
};
