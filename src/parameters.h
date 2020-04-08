#include <utility>

//
// Created by krab1k on 05/11/18.
//

#pragma once

#include <map>
#include <memory>
#include <utility>
#include <vector>
#include <functional>
#include <tuple>
#include <string>

#include "structures/atom.h"
#include "structures/bond.h"
#include "method.h"


class MoleculeSet;


class CommonParameters {
    friend class Parameters;

    std::vector<std::string> names_;
    std::vector<double> parameters_;
public:
    explicit CommonParameters(std::vector<std::string> names, std::vector<double> parameters) :
            names_{std::move(names)}, parameters_{std::move(parameters)} {}

    [[nodiscard]] const std::vector<std::string> &names() const { return names_; }

    [[nodiscard]] double parameter(size_t idx) const { return parameters_[idx]; }
};

class AtomParameters {
    friend class Parameters;

    std::vector<std::string> names_;
    std::vector<std::tuple<std::string, std::string, std::string>> keys_;
    std::vector<std::vector<double>> parameters_;
public:
    explicit AtomParameters(std::vector<std::string> names,
                            std::vector<std::vector<double>> parameters,
                            std::vector<std::tuple<std::string, std::string, std::string>> parameter_order)
            : names_{std::move(names)}, keys_{std::move(parameter_order)},
              parameters_{std::move(parameters)} {}

    [[nodiscard]] const std::vector<std::string> &names() const { return names_; }

    [[nodiscard]] const std::vector<std::tuple<std::string, std::string, std::string>> &keys() const { return keys_; }

    [[nodiscard]] std::function<double(const Atom &)> parameter(size_t idx) const noexcept;
};

class BondParameters {
    friend class Parameters;

    std::vector<std::string> names_;
    std::vector<std::tuple<std::string, std::string, std::string, std::string>> keys_;
    std::vector<std::vector<double>> parameters_;
public:
    explicit BondParameters(std::vector<std::string> names,
                            std::vector<std::vector<double>> parameters,
                            std::vector<std::tuple<std::string, std::string, std::string, std::string>> parameter_order)
            : names_{std::move(names)}, keys_{std::move(parameter_order)},
              parameters_{std::move(parameters)} {}

    [[nodiscard]] const std::vector<std::string> &names() const { return names_; }

    [[nodiscard]] const std::vector<std::tuple<std::string, std::string, std::string, std::string>> &keys() const { return keys_; }

    [[nodiscard]] std::function<double(const Bond &)> parameter(size_t idx) const noexcept;
};

class Parameters {
    std::string name_{};
    std::string method_name_{};
    std::string source_{};

    std::unique_ptr<CommonParameters> common_{nullptr};
    std::unique_ptr<AtomParameters> atoms_{nullptr};
    std::unique_ptr<BondParameters> bonds_{nullptr};

public:
    explicit Parameters(const std::string &filename);

    explicit Parameters(const MoleculeSet &ms, const Method *method);

    void print() const;

    void save_to_file(const std::string &filename) const;

    [[nodiscard]] const std::string &name() const { return name_; }

    [[nodiscard]] const std::string &method_name() const { return method_name_; }

    [[nodiscard]] const std::string &source() const { return source_; }

    [[nodiscard]] std::vector<double> get_vector() const;

    void set_from_vector(const std::vector<double> &parameters);

    [[nodiscard]] const CommonParameters *common() const { return common_.get(); }

    [[nodiscard]] const AtomParameters *atom() const { return atoms_.get(); }

    [[nodiscard]] const BondParameters *bond() const { return bonds_.get(); }
};
