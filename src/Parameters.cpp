//
// Created by krab1k on 05/11/18.
//

#include "Parameters.h"

#include <iostream>
#include <map>
#include <vector>
#include <QString>
#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>


Parameters::Parameters(const QString &filename) {
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        std::cerr << "Cannot open file: " << filename.toStdString() << std::endl;
        exit(EXIT_FAILURE);
    }

    QByteArray json_data = file.readAll();
    QJsonObject document = QJsonDocument::fromJson(json_data).object();
    if (document.contains("common")) {
        QJsonObject common = document["common"].toObject();
        QJsonArray names_json = common["names"].toArray();
        QJsonArray values_json = common["values"].toArray();
        if (names_json.size() != values_json.size()) {
            std::cerr << "Incorrect parameter file" << std::endl;
            exit(EXIT_FAILURE);
        }

        std::vector<QString> names;

        std::map<QString, double> parameters;
        for (int i = 0; i < names_json.size(); i++) {
            names.emplace_back(names_json[i].toString());
            parameters[names_json[i].toString()] = values_json[i].toDouble();
        }

        common_ = std::make_unique<CommonParameters>(names, parameters);
    }

    if (document.contains("atom")) {
        QJsonObject atom = document["atom"].toObject();
        QJsonArray names_json = atom["names"].toArray();
        QJsonArray data_json = atom["data"].toArray();

        std::vector<QString> names;
        for(auto &&name: names_json) {
            names.emplace_back(name.toString());
        }

        std::vector<std::tuple<QString, QString, QString> > parameter_order;
        std::map<std::tuple<QString, QString, QString>, std::vector<double>> parameters;
        for (auto &&ref_pair: data_json) {
            QJsonObject pair = ref_pair.toObject();
            QJsonArray json_key = pair["key"].toArray();
            QJsonArray json_values = pair["value"].toArray();

            if (names_json.size() != json_values.size()) {
                std::cerr << "Incorrect parameter file" << std::endl;
                exit(EXIT_FAILURE);
            }

            auto key = std::make_tuple(json_key[0].toString(), json_key[1].toString(), json_key[2].toString());
            parameter_order.push_back(key);

            std::vector<double> values;
            for (auto &&value: json_values) {
                values.push_back(value.toDouble());
            }
            parameters[key] = values;
        }

        atoms_ = std::make_unique<AtomParameters>(names, parameters, parameter_order);
    }
    if (document.contains("bond")) {
        QJsonObject bond = document["bond"].toObject();
        QJsonArray names_json = bond["names"].toArray();
        QJsonArray data_json = bond["data"].toArray();

        std::vector<QString> names;
        std::vector<std::tuple<QString, QString, QString, QString>> parameter_order;
        std::map<std::tuple<QString, QString, QString, QString>, std::vector<double>> parameters;

        for(auto &&name: names_json) {
            names.emplace_back(name.toString());
        }

        for (auto &&ref_pair: data_json) {
            QJsonObject pair = ref_pair.toObject();
            QJsonArray json_key = pair["key"].toArray();
            QJsonArray json_values = pair["value"].toArray();

            if (names_json.size() != json_values.size()) {
                std::cerr << "Incorrect parameter file" << std::endl;
                exit(EXIT_FAILURE);
            }

            auto key = std::make_tuple(json_key[0].toString(), json_key[1].toString(), json_key[2].toString(),
                                       json_key[3].toString());

            parameter_order.push_back(key);

            std::vector<double> values;
            for (auto &&value: json_values) {
                values.push_back(value.toDouble());
            }
            parameters[key] = values;
        }

        bonds_ = std::make_unique<BondParameters>(names, parameters, parameter_order);
    }
}

void Parameters::print() const {
    if (common_) {
        std::cout << "Common parameters" << std::endl;
        for (auto &[key, val]: common_->parameters_) {
            std::cout << key.toStdString() << ": " << val << std::endl;
        }
    }
    if (atoms_) {
        std::cout << "Atom parameters" << std::endl;
        for (const auto &key: atoms_->parameter_order_) {
            auto &[symbol, cls, type] = key;
            std::cout << symbol.toStdString() << " " << cls.toStdString() << " " << type.toStdString() << ": ";
            for (double val: atoms_->parameters_[key]) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }
    if (bonds_) {
        std::cout << "Bond parameters" << std::endl;
        for (const auto &key: bonds_->parameter_order_) {
            auto &[symbol1, symbol2, cls, type] = key;
            std::cout << symbol1.toStdString() << " " << symbol2.toStdString() << " " << cls.toStdString() << " "
                      << type.toStdString() << ": ";

            for (double val: bonds_->parameters_[key]) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }
}

std::function<double(const Atom &)> AtomParameters::parameter(const QString &name) const {
    long idx;
    auto it = std::find(names_.begin(), names_.end(), name);
    if (it != names_.end()) {
        idx = std::distance(names_.begin(), it);
    } else {
        std::cerr << "Invalid parameter name: " << name.toStdString() << std::endl;
        exit(EXIT_FAILURE);
    }
    return [this, idx](const Atom &atom) { return parameters_.at(atom.atom_type())[idx]; };
}