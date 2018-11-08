//
// Created by krab1k on 05/11/18.
//

#include "Parameters.h"

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

Parameters::Parameters(const std::string &filename) {
    namespace pt = boost::property_tree;

    pt::ptree iroot;
    pt::read_json(filename, iroot);

    auto it = iroot.find("common");
    if (it != iroot.not_found()){
        std::vector<std::string> names;
        std::vector<double> values;
        for(const auto &name: iroot.get_child("common.names")) {
            names.push_back(name.second.data());
        }
        for(const auto &value: iroot.get_child("common.values")) {
            values.push_back(value.second.get_value<double>());
        }
        common_ = std::make_unique<CommonParameters>(names, values);
    }

    it = iroot.find("atom");
    if (it != iroot.not_found()) {
        std::vector<std::string> names;
        std::vector<std::tuple<std::string, std::string, std::string>> keys;
        std::vector<std::vector<double>> parameters;

        for(const auto &name: iroot.get_child("atom.names")) {
            names.push_back(name.second.data());
        }

        for(const auto &pair: iroot.get_child("atom.data")) {
            std::vector<std::string> tmp_types;
            std::vector<double> values;
            for(const auto &t: pair.second.get_child("key")) {
                tmp_types.push_back(t.second.data());
            }
            keys.emplace_back(std::make_tuple(tmp_types[0], tmp_types[1], tmp_types[2]));
            tmp_types.clear();
            for(const auto &v: pair.second.get_child("value")) {
                values.push_back(v.second.get_value<double>());
            }
            parameters.push_back(values);
        }
        atoms_ = std::make_unique<AtomParameters>(names, parameters, keys);
    }

    it = iroot.find("bond");
    if (it != iroot.not_found()) {
        std::vector<std::string> names;
        std::vector<std::tuple<std::string, std::string, std::string, std::string>> keys;
        std::vector<std::vector<double>> parameters;

        for(const auto &name: iroot.get_child("bond.names")) {
            names.push_back(name.second.data());
        }

        for(const auto &pair: iroot.get_child("bond.data")) {
            std::vector<std::string> tmp_types;
            std::vector<double> values;
            for(const auto &t: pair.second.get_child("key")) {
                tmp_types.push_back(t.second.data());
            }
            keys.emplace_back(std::make_tuple(tmp_types[0], tmp_types[1], tmp_types[2], tmp_types[3]));
            tmp_types.clear();
            for(const auto &v: pair.second.get_child("value")) {
                values.push_back(v.second.get_value<double>());
            }
            parameters.push_back(values);
        }
        bonds_ = std::make_unique<BondParameters>(names, parameters, keys);
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