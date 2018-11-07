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

#include "structures/Atom.h"

class CommonParameters {
    friend class Parameters;

    std::vector<std::string> names_;
    std::map<std::string, double> parameters_;
public:
    explicit CommonParameters(std::vector<std::string> names, std::map<std::string, double> parameters) :
            names_{std::move(names)}, parameters_{std::move(parameters)} {}

    const std::vector<std::string> &names() const { return names_; }

    double parameter(const std::string &name) const { return parameters_.at(name); }

};

class AtomParameters {
    friend class Parameters;

    std::vector<std::string> names_;
    std::vector<std::tuple<std::string, std::string, std::string>> keys_;
    std::map<std::tuple<std::string, std::string, std::string>, std::vector<double>> parameters_;
public:
    explicit AtomParameters(std::vector<std::string> names,
                            std::map<std::tuple<std::string, std::string, std::string>, std::vector<double> > parameters,
                            std::vector<std::tuple<std::string, std::string, std::string>> parameter_order)
            : names_{std::move(names)}, keys_{std::move(parameter_order)},
              parameters_{std::move(parameters)} {}

    const std::vector<std::string> &names() const { return names_; }

    const std::vector<std::tuple<std::string, std::string, std::string>> &keys() const { return keys_; }

    std::function<double(const Atom &)> parameter(const std::string &name) const;


};

class BondParameters {
    friend class Parameters;

    std::vector<std::string> names_;
    std::vector<std::tuple<std::string, std::string, std::string, std::string>> keys_;
    std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<double>> parameters_;
public:
    explicit BondParameters(std::vector<std::string> names,
                            std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<double>> parameters,
                            std::vector<std::tuple<std::string, std::string, std::string, std::string>> parameter_order)
            : names_{std::move(names)}, keys_{std::move(parameter_order)},
              parameters_{std::move(parameters)} {}

    const std::vector<std::string> &names() const { return names_; }

    const std::vector<std::tuple<std::string, std::string, std::string, std::string>> &keys() const { return keys_; }

};

class Parameters {
    std::unique_ptr<CommonParameters> common_{nullptr};
    std::unique_ptr<AtomParameters> atoms_{nullptr};
    std::unique_ptr<BondParameters> bonds_{nullptr};

public:
    explicit Parameters(const std::string &filename);

    void print() const;

    const CommonParameters *common() const { return common_.get(); }

    const AtomParameters *atom() const { return atoms_.get(); }

    const BondParameters *bond() const { return bonds_.get(); }
};






