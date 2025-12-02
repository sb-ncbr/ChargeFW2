#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>
#include <utility>

#include "structures/molecule.h"


enum class RequiredFeatures {
    BOND_INFO,
    BOND_DISTANCES,
    DISTANCE_TREE
};


class Parameters;


struct MethodOption {
    std::string name;
    std::string description;
    std::string type;
    std::string default_value;
    std::vector<std::string> choices;
};

struct MethodMetadata {
    std::string name;
    std::string internal_name;
    std::string full_name;
    std::optional<std::string> publication;
    std::string type;
    uint16_t priority;
};

class Method {
protected:
    const std::vector<std::string> common_parameters_{};
    const std::vector<std::string> atom_parameters_{};
    const std::vector<std::string> bond_parameters_{};
    const std::map<std::string, MethodOption> options_{};

    std::map<std::string, std::string> option_values_{};

    Parameters *parameters_{nullptr};

public:
    Method(std::vector<std::string> common, std::vector<std::string> atom,
           std::vector<std::string> bond, std::map<std::string, MethodOption> options) :
            common_parameters_{std::move(common)},
            atom_parameters_{std::move(atom)},
            bond_parameters_{std::move(bond)},
            options_{std::move(options)} {}

    virtual ~Method() = default;

    [[nodiscard]] const std::vector<std::string> &common_parameters() const { return common_parameters_; }

    [[nodiscard]] const std::vector<std::string> &atom_parameters() const { return atom_parameters_; }

    [[nodiscard]] const std::vector<std::string> &bond_parameters() const { return bond_parameters_; }

    Parameters *parameters() { return parameters_; }

    void set_parameters(Parameters *parameters);

    [[nodiscard]] bool has_parameters() const {
        return (common_parameters_.size() + atom_parameters_.size() + bond_parameters_.size()) != 0;
    }

    [[nodiscard]] virtual std::vector<RequiredFeatures> get_requirements() const;

    [[nodiscard]] virtual bool is_suitable_for_molecule(const Molecule &) const;

    [[nodiscard]] virtual bool is_suitable_for_large_molecule() const { return true; }

    [[nodiscard]] virtual std::vector<double> calculate_charges(const Molecule &molecule) const = 0;

    [[nodiscard]] std::map<std::string, MethodOption> get_options() const { return options_; }

    [[nodiscard]] virtual const MethodMetadata& metadata() const = 0;

    template<typename T>
    T get_option_value(const std::string &name) const;

    void set_option_value(const std::string &name, const std::string &value) {
        option_values_[name] = value;
    }
};


template<>
[[nodiscard]] std::string Method::get_option_value<std::string>(const std::string &name) const;

template<>
[[nodiscard]] double Method::get_option_value<double>(const std::string &name) const;

template<>
[[nodiscard]] int Method::get_option_value<int>(const std::string &name) const;

std::unique_ptr<Method> load_method(std::string const& name);

std::vector<std::unique_ptr<Method>> get_available_methods();

template<class T>
std::unique_ptr<Method> make_method() { return std::make_unique<T>(); }
