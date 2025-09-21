#include <memory>
#include <regex>
#include <print>
#include <format>
#include <vector>
#include <map>
#include <set>
#include <filesystem>
#include <dlfcn.h>
#include <Eigen/Dense>
#include <omp.h>

#include "method.h"
#include "parameters.h"
#include "exceptions/file_exception.h"
#include "utility/install.h"


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


Method* load_method(const std::string &method_name) {
    std::string file = InstallPaths::libdir() / ("lib" + method_name + ".so");

    auto handle = dlopen(file.c_str(), RTLD_LAZY);

    auto get_method_handle = reinterpret_cast<Method *(*)()>(dlsym(handle, "get_method"));
    if (!get_method_handle) {
        throw FileException(dlerror());
    }

    return (*get_method_handle)();
}


std::vector<Method*> get_available_methods() {
    std::vector<Method*> results;
    std::regex method_pattern(R"(^lib(.*)\.so$)");

    for (const auto &entry : fs::directory_iterator(InstallPaths::libdir())) {
        auto filename = entry.path().filename().string();
        std::smatch matches;

        if (std::regex_match(filename, matches, method_pattern)) {
            auto method_name = matches[1].str();
            Method* method;
            
            try {
                method = load_method(method_name);
            } catch (FileException &e){
                std::println(stderr, "Failed to load method: {}", method_name);
                std::println("{}", e.what());
                continue;
            }
            
            results.emplace_back(method);
        }
    }

    std::ranges::sort(results, [](const auto &a, const auto &b) {
        return a->metadata().priority > b->metadata().priority;
    });

    return results;
}
