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
    explicit Charges(std::string method_name = {}, std::string parameters_name = "None")
        : method_name_(std::move(method_name)), parameters_name_(std::move(parameters_name)) {
    }

    [[nodiscard]] std::string method_name() const { return method_name_; }

    [[nodiscard]] std::string parameters_name() const { return parameters_name_; }

    [[nodiscard]] std::vector<std::string> names() const { return names_; }

    std::vector<double> operator[](const std::string& name) const { return charges_.at(name); }

    void insert(const std::string& name, std::vector<double> charges);
};
