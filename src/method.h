//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <string>
#include <vector>

#include "structures/molecule.h"


class Parameters;

struct MethodOption {
    std::string name;
    std::string description;
    std::string type;
    std::string default_value;
    std::vector<std::string> choices;
};

class Method {
protected:
    const std::string name_{};
    const std::vector<std::string> common_parameters_{};
    const std::vector<std::string> atom_parameters_{};
    const std::vector<std::string> bond_parameters_{};
    const std::map<std::string, MethodOption> options_{};

    std::map<std::string, std::string> option_values_{};

    Parameters *parameters_{nullptr};

public:
    Method(std::string name, std::vector<std::string> common, std::vector<std::string> atom,
           std::vector<std::string> bond, std::map<std::string, MethodOption> options) :
            name_{std::move(name)},
            common_parameters_{std::move(common)},
            atom_parameters_{std::move(atom)},
            bond_parameters_{std::move(bond)},
            options_{std::move(options)} {}

    std::vector<std::string> common_parameters() { return common_parameters_; }

    std::vector<std::string> atom_parameters() { return atom_parameters_; }

    std::vector<std::string> bond_parameters() { return bond_parameters_; }

    Parameters *parameters() { return parameters_; }

    void set_parameters(Parameters *parameters);

    bool has_parameters() {
        return (common_parameters_.size() + atom_parameters_.size() + bond_parameters_.size()) != 0;
    }

    virtual std::vector<double> calculate_charges(const Molecule &molecule) const = 0;

    std::string name() { return name_; }

    std::map<std::string, MethodOption> get_options() const { return options_; }

    template<typename T>
    T get_option_value(const std::string &name) const;

    void set_option_value(const std::string &name, const std::string &value) {
        option_values_[name] = value;
    }
};

