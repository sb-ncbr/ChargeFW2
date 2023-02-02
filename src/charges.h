//
// Created by krab1k on 8.11.18.
//

#pragma once

#include <string>
#include <vector>
#include <map>


class Charges {
    std::string method_name_{};
    std::string parameters_name_{"None"};
    std::vector<std::string> names_{};
    std::map<std::string, std::vector<double>> charges_{};
public:
    Charges() = default;

    explicit Charges(const std::string &filename);

    explicit Charges(const std::map<std::string, std::vector<double>> &charges);

    void set_method_name(const std::string &method_name) { method_name_ = method_name; }

    void set_parameters_name(const std::string &parameters_name) {parameters_name_ = parameters_name; }

    [[nodiscard]] std::string method_name() const { return method_name_; }

    [[nodiscard]] std::string parameters_name() const { return parameters_name_; }

    [[nodiscard]] std::vector<std::string> names() const { return names_; }

    std::vector<double> operator[](const std::string &name) const { return charges_.at(name); }

    void insert(const std::string &name, std::vector<double> charges);
};
