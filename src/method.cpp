//
// Created by krab1k on 6.11.18.
//

#include <fmt/format.h>
#include <fmt/ranges.h>

#include "method.h"
#include "parameters.h"
#include "config.h"


void Method::set_parameters(Parameters *parameters) {
    if (common_parameters_.size() + atom_parameters_.size() + bond_parameters_.size() == 0 and parameters == nullptr) {
        return;
    }
    if (parameters->common() != nullptr and parameters->common()->names() != common_parameters_) {
        fmt::print(stderr, "Invalid common parameters provided\n");
        fmt::print(stderr, "Expected: {}\n", common_parameters_);
        fmt::print(stderr, "Got: {}\n", parameters->common()->names());
        exit(EXIT_FILE_ERROR);
    }

    if (parameters->atom() != nullptr and parameters->atom()->names() != atom_parameters_) {
        fmt::print(stderr, "Invalid atom parameters provided\n");
        fmt::print(stderr, "Expected: {}\n", atom_parameters_);
        fmt::print(stderr, "Got: {}\n", parameters->atom()->names());
        exit(EXIT_FILE_ERROR);
    }

    if (parameters->bond() != nullptr and parameters->bond()->names() != bond_parameters_) {
        fmt::print(stderr, "Invalid bond parameters provided\n");
        fmt::print(stderr, "Expected: {}\n", bond_parameters_);
        fmt::print(stderr, "Got: {}\n", parameters->bond()->names());
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


std::vector<double> EEMethod::calculate_charges(const Molecule &molecule) const {

    auto method = get_option_value<std::string>("type");
    auto radius = get_option_value<double>("radius");
    std::vector<const Atom *> fragment_atoms;
    if (method == "full") {
        for (const auto &atom: molecule.atoms()) {
            fragment_atoms.push_back(&atom);
        }

        return solve_system(fragment_atoms, molecule.total_charge());

    } else /* method == "cutoff" */ {
        std::vector<double> results;
        for (const auto &atom: molecule.atoms()) {
            fragment_atoms = molecule.get_close_atoms(atom, radius);
            auto res = solve_system(fragment_atoms,
                                        static_cast<double>(molecule.total_charge()) * fragment_atoms.size() /
                                        molecule.atoms().size());
            results.push_back(res[0]);
        }

        double correction = molecule.total_charge();
        for (auto val: results) {
            correction -= val;
        }

        correction /= molecule.atoms().size();

        for (auto &val: results) {
            val += correction;
        }

        return results;
    }
}
