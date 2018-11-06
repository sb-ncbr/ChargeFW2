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
#include<QString>

#include "structures/Atom.h"

class CommonParameters {
    friend class Parameters;

    std::vector<QString> names_;
    std::map<QString, double> parameters_;
public:
    explicit CommonParameters(std::vector<QString> names, std::map<QString, double> parameters) :
            names_{std::move(names)}, parameters_{std::move(parameters)} {}

    const std::vector<QString> &names() const { return names_; }

    double parameter(const QString &name) const { return parameters_.at(name); }

};

class AtomParameters {
    friend class Parameters;

    std::vector<QString> names_;
    std::vector<std::tuple<QString, QString, QString>> parameter_order_;
    std::map<std::tuple<QString, QString, QString>, std::vector<double>> parameters_;
public:
    explicit AtomParameters(std::vector<QString> names,
                            std::map<std::tuple<QString, QString, QString>, std::vector<double> > parameters,
                            std::vector<std::tuple<QString, QString, QString>> parameter_order)
            : names_{std::move(names)}, parameter_order_{std::move(parameter_order)},
              parameters_{std::move(parameters)} {}

    const std::vector<QString> &names() const { return names_; }

    const std::vector<std::tuple<QString, QString, QString>> &keys() const { return parameter_order_; }

    std::function<double(const Atom &)> parameter(const QString &name) const;


};

class BondParameters {
    friend class Parameters;

    std::vector<QString> names_;
    std::vector<std::tuple<QString, QString, QString, QString>> parameter_order_;
    std::map<std::tuple<QString, QString, QString, QString>, std::vector<double>> parameters_;
public:
    explicit BondParameters(std::vector<QString> names,
                            std::map<std::tuple<QString, QString, QString, QString>, std::vector<double>> parameters,
                            std::vector<std::tuple<QString, QString, QString, QString>> parameter_order)
            : names_{std::move(names)}, parameter_order_{std::move(parameter_order)},
              parameters_{std::move(parameters)} {}

    const std::vector<QString> &names() const { return names_; }

    const std::vector<std::tuple<QString, QString, QString, QString>> &keys() const { return parameter_order_; }

};

class Parameters {
    std::unique_ptr<CommonParameters> common_{nullptr};
    std::unique_ptr<AtomParameters> atoms_{nullptr};
    std::unique_ptr<BondParameters> bonds_{nullptr};

public:
    explicit Parameters(const QString &filename);

    void print() const;

    const CommonParameters *common() const { return common_.get(); }

    const AtomParameters *atom() const { return atoms_.get(); }

    const BondParameters *bond() const { return bonds_.get(); }
};






