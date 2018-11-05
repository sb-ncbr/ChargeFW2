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
    std::map<QString, double> parameters_;
public:
    explicit CommonParameters(std::map<QString, double> parameters) : parameters_{std::move(parameters)} {}

    double operator[](const QString &name) { return parameters_.at(name); }

    friend class Parameters;
};

class AtomParameters {
    std::vector<QString> names_;
    std::vector<std::tuple<QString, QString, QString>> parameter_order_;
    std::map<std::tuple<QString, QString, QString>, std::vector<double>> parameters_;
public:
    explicit AtomParameters(std::vector<QString> names,
                            std::map<std::tuple<QString, QString, QString>, std::vector<double> > parameters,
                            std::vector<std::tuple<QString, QString, QString>> parameter_order)
            : names_{std::move(names)}, parameter_order_{std::move(parameter_order)},
              parameters_{std::move(parameters)} {}

    const std::vector<std::tuple<QString, QString, QString>> &keys() const { return parameter_order_; }

    std::function<double(const Atom &)> operator[](const QString &name) const;

    friend class Parameters;

};

class BondParameters {
    std::vector<QString> names_;
    std::vector<std::tuple<QString, QString, QString, QString>> parameter_order_;
    std::map<std::tuple<QString, QString, QString, QString>, std::vector<double>> parameters_;
public:
    explicit BondParameters(std::vector<QString> names,
                            std::map<std::tuple<QString, QString, QString, QString>, std::vector<double>> parameters,
                            std::vector<std::tuple<QString, QString, QString, QString>> parameter_order)
            : names_{std::move(names)}, parameter_order_{std::move(parameter_order)},
              parameters_{std::move(parameters)} {}

    const std::vector<std::tuple<QString, QString, QString, QString>> &keys() const { return parameter_order_; }

    friend class Parameters;
};

class Parameters {
    std::unique_ptr<CommonParameters> common_{nullptr};
    std::unique_ptr<AtomParameters> atoms_{nullptr};
    std::unique_ptr<BondParameters> bonds_{nullptr};

public:
    explicit Parameters(const QString &filename);

    void print() const;

    const CommonParameters &common() const { return *common_; }

    const AtomParameters &atom() const { return *atoms_; }

    const BondParameters &bond() const { return *bonds_; }
};






