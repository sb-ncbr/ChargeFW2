#include <memory>
#include <regex>
#include <print>
#include <format>
#include <vector>
#include <map>
#include <algorithm>
#include <filesystem>

#include "method.h"
#include "parameters.h"
#include "utility/install.h"
#include "utility/exceptions.h"
#include "method_registry.h"


namespace fs = std::filesystem;


std::vector<RequiredFeatures> Method::get_requirements() const {
    return {};
}


void Method::set_parameters(Parameters *parameters) {
    if (common_parameters_.size() + atom_parameters_.size() + bond_parameters_.size() == 0 and parameters == nullptr) {
        return;
    }
    if (parameters->common() != nullptr and parameters->common()->names() != common_parameters_) {
        std::println(stderr, "Parameters don't match");
        std::println(stderr, "Expected: {}\n", common_parameters_);
        std::println(stderr, "Got: {}\n", parameters->common()->names());
        throw std::runtime_error("Invalid common parameters provided");
    }

    if (parameters->atom() != nullptr and parameters->atom()->names() != atom_parameters_) {
        std::println(stderr, "Parameters don't match");
        std::println(stderr, "Expected: {}", atom_parameters_);
        std::println(stderr, "Got: {}", parameters->atom()->names());
        throw std::runtime_error("Invalid atom parameters provided");
    }

    if (parameters->bond() != nullptr and parameters->bond()->names() != bond_parameters_) {
        std::println(stderr, "Parameters don't match");
        std::println(stderr, "Expected: {}", bond_parameters_);
        std::println(stderr, "Got: {}", parameters->bond()->names());
        throw std::runtime_error("Invalid bond parameters provided");
    }
    parameters_ = parameters;
}


template<>
std::string Method::get_option_value<std::string>(const std::string &name) const {
    return option_values_.at(name);
}


bool Method::is_suitable_for_molecule(const Molecule &) const {
    return true;
}


template<>
double Method::get_option_value<double>(const std::string &name) const {
    return std::stod(option_values_.at(name));
}


template<>
int Method::get_option_value<int>(const std::string &name) const {
    return std::stoi(option_values_.at(name));
}


std::unique_ptr<Method> load_method(std::string const& name) {
    auto m = MethodRegistry::make(name);     // already returns unique_ptr<Method>
    if (!m) throw FileException("Unknown method: " + name);
    return m;                                // no release() — ownership stays with caller safely
}


std::vector<std::unique_ptr<Method>> get_available_methods() {
    std::vector<std::unique_ptr<Method>> v;
    v.reserve(MethodRegistry::names().size());

    for (auto const& name : MethodRegistry::names()) {
        if (auto m = MethodRegistry::make(name)) v.emplace_back(std::move(m));
    }

    std::ranges::sort(v, [](auto const& a, auto const& b) {
        return a->metadata().priority > b->metadata().priority;
    });

    return v;
}
