//
// Created by krab1k on 6.11.18.
//

#include <iostream>

#include "method.h"
#include "parameters.h"
#include "utility/utility.h"
#include "config.h"


void Method::set_parameters(const Parameters *parameters) {
    if (common_parameters_.size() + atom_parameters_.size() + bond_parameters_.size() == 0 and parameters == nullptr) {
        return;
    }
    if (parameters->common() != nullptr and parameters->common()->names() != common_parameters_) {
        std::cerr << "Invalid common parameters provided" << std::endl;
        std::cerr << "Expected: " << common_parameters_ << std::endl;
        std::cerr << "Got: " << parameters->common()->names() << std::endl;
        exit(EXIT_FILE_ERROR);
    }

    if (parameters->atom() != nullptr and parameters->atom()->names() != atom_parameters_) {
        std::cerr << "Invalid atom parameters provided" << std::endl;
        std::cerr << "Expected: " << atom_parameters_ << std::endl;
        std::cerr << "Got: " << parameters->atom()->names() << std::endl;
        exit(EXIT_FILE_ERROR);
    }

    if (parameters->bond() != nullptr and parameters->bond()->names() != bond_parameters_) {
        std::cerr << "Invalid bond parameters provided" << std::endl;
        std::cerr << "Expected: " << bond_parameters_ << std::endl;
        std::cerr << "Got: " << parameters->bond()->names() << std::endl;
        exit(EXIT_FILE_ERROR);
    }
    parameters_ = parameters;
}


template<>
std::string Method::get_option_value<std::string>(const std::string &name) const {
    return option_values_.at(name);
}


template<>
double Method::get_option_value<double>(const std::string &name) const {
    return std::stod(option_values_.at(name));
}


template<>
int Method::get_option_value<int>(const std::string &name) const {
    return std::stoi(option_values_.at(name));
}