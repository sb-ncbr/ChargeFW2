//
// Created by krab1k on 05/11/18.
//

#include "Parameters.h"

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <QString>
#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>


Parameters::Parameters(const std::string &filename) {
    QFile file(QString::fromStdString(filename));
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        std::cerr << "Cannot open file: " << filename << std::endl;
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

        std::vector<std::string> names;

        std::vector<double> parameters;
        for (int i = 0; i < names_json.size(); i++) {
            names.emplace_back(names_json[i].toString().toStdString());
            parameters.push_back(values_json[i].toDouble());
        }

        common_ = std::make_unique<CommonParameters>(names, parameters);
    }

    if (document.contains("atom")) {
        QJsonObject atom = document["atom"].toObject();
        QJsonArray names_json = atom["names"].toArray();
        QJsonArray data_json = atom["data"].toArray();

        std::vector<std::string> names;
        for (auto &&name: names_json) {
            names.emplace_back(name.toString().toStdString());
        }

        std::vector<std::tuple<std::string, std::string, std::string>> parameter_order;
        std::vector<std::vector<double>> parameters;
        for (auto &&ref_pair: data_json) {
            QJsonObject pair = ref_pair.toObject();
            QJsonArray json_key = pair["key"].toArray();
            QJsonArray json_values = pair["value"].toArray();

            if (names_json.size() != json_values.size()) {
                std::cerr << "Incorrect parameter file" << std::endl;
                exit(EXIT_FAILURE);
            }

            auto key = std::make_tuple(json_key[0].toString().toStdString(), json_key[1].toString().toStdString(),
                                       json_key[2].toString().toStdString());
            parameter_order.push_back(key);

            std::vector<double> values;
            for (auto &&value: json_values) {
                values.push_back(value.toDouble());
            }
            parameters.push_back(values);
        }

        atoms_ = std::make_unique<AtomParameters>(names, parameters, parameter_order);
    }
    if (document.contains("bond")) {
        QJsonObject bond = document["bond"].toObject();
        QJsonArray names_json = bond["names"].toArray();
        QJsonArray data_json = bond["data"].toArray();

        std::vector<std::string> names;
        std::vector<std::tuple<std::string, std::string, std::string, std::string>> parameter_order;
        std::vector<std::vector<double>> parameters;

        for (auto &&name: names_json) {
            names.emplace_back(name.toString().toStdString());
        }

        for (auto &&ref_pair: data_json) {
            QJsonObject pair = ref_pair.toObject();
            QJsonArray json_key = pair["key"].toArray();
            QJsonArray json_values = pair["value"].toArray();

            if (names_json.size() != json_values.size()) {
                std::cerr << "Incorrect parameter file" << std::endl;
                exit(EXIT_FAILURE);
            }

            auto key = std::make_tuple(json_key[0].toString().toStdString(), json_key[1].toString().toStdString(),
                                       json_key[2].toString().toStdString(), json_key[3].toString().toStdString());

            parameter_order.push_back(key);

            std::vector<double> values;
            for (auto &&value: json_values) {
                values.push_back(value.toDouble());
            }
            parameters.push_back(values);
        }

        bonds_ = std::make_unique<BondParameters>(names, parameters, parameter_order);
    }
}

void Parameters::print() const {
    if (common_) {
        std::cout << "Common parameters" << std::endl;
        for(size_t i = 0; i < common_->names().size(); i++) {
            std::cout << common_->names()[i] << ": " << common_->parameters_[i] << std::endl;
        }
    }
    if (atoms_) {
        std::cout << "Atom parameters" << std::endl;
        for (size_t i = 0; i < atoms_->parameters_.size(); i++) {
            auto &[symbol, cls, type] = atoms_->keys()[i];
            std::cout << symbol << " " << cls << " " << type << ": ";
            for (double val: atoms_->parameters_[i]) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }
    if (bonds_) {
        std::cout << "Bond parameters" << std::endl;
        for (size_t i = 0; i < bonds_->parameters_.size(); i++) {
            auto &[symbol1, symbol2, cls, type] = bonds_->keys()[i];
            std::cout << symbol1 << " " << symbol2 << " " << cls << " " << type << ": ";

            for (double val: bonds_->parameters_[i]) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }
}

std::function<double(const Atom &)> AtomParameters::parameter(int idx) const {

    return [this, idx](const Atom &atom) { return parameters_[atom.atom_type()][idx]; };
}