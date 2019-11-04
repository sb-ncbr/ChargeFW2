//
// Created by krab1k on 05/11/18.
//

#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <nlohmann/json.hpp>
#include <fmt/format.h>

#include "chargefw2.h"
#include "parameters.h"
#include "structures/molecule_set.h"
#include "method.h"


Parameters::Parameters(const std::string &filename) {
    using json = nlohmann::json;
    json j;
    std::ifstream f(filename);
    if (!f) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    f >> j;
    f.close();

    try {
        name_ = j["metadata"]["name"];
        method_name_ = j["metadata"]["method"];

        if (j["metadata"].count("source")) {
            source_ = j["metadata"]["source"];
        } else {
            source_ = "small";
        }

        if (j.count("common")) {
            auto names = j["common"]["names"].get<std::vector<std::string>>();
            auto values = j["common"]["values"].get<std::vector<double>>();
            common_ = std::make_unique<CommonParameters>(names, values);
        }

        if (j.count("atom")) {
            auto names = j["atom"]["names"].get<std::vector<std::string>>();
            std::vector<std::tuple<std::string, std::string, std::string>> keys;
            std::vector<std::vector<double>> parameters;
            for (const auto &obj: j["atom"]["data"]) {
                auto key = obj["key"];
                keys.emplace_back(
                        std::make_tuple(key[0].get<std::string>(), key[1].get<std::string>(), key[2].get<std::string>()));
                parameters.emplace_back(obj["value"].get<std::vector<double>>());
            }
            atoms_ = std::make_unique<AtomParameters>(names, parameters, keys);
        }

        if (j.count("bond")) {
            auto names = j["bond"]["names"].get<std::vector<std::string>>();
            std::vector<std::tuple<std::string, std::string, std::string, std::string>> keys;
            std::vector<std::vector<double>> parameters;
            for (const auto &obj: j["bond"]["data"]) {
                auto key = obj["key"];
                keys.emplace_back(
                        std::make_tuple(key[0].get<std::string>(), key[1].get<std::string>(), key[2].get<std::string>(),
                                        key[3].get<std::string>()));
                parameters.emplace_back(obj["value"].get<std::vector<double>>());
            }
            bonds_ = std::make_unique<BondParameters>(names, parameters, keys);
        }
    }
    catch (std::exception &) {
        fmt::print(stderr, "Incorrect file with parameters\n");
        exit(EXIT_FILE_ERROR);
    }
}


void Parameters::save_to_file(const std::string &filename) const {
    using json = nlohmann::json;

    json j;
    j["metadata"]["name"] = name_;
    j["metadata"]["method"] = method_name_;
    j["metadata"]["publication"] = nullptr;

    if (common_) {
        j["common"]["names"] = common_->names_;
        j["common"]["values"] = common_->parameters_;
    }

    if (atoms_) {
        j["atom"]["names"] = atoms_->names_;
        for (size_t i = 0; i < atoms_->keys_.size(); i++) {
            auto obj = json::object();
            obj["key"] = atoms_->keys_[i];
            obj["value"] = atoms_->parameters_[i];
            j["atom"]["data"].push_back(obj);
        }
    }

    if (bonds_) {
        j["bond"]["names"] = bonds_->names_;
        for (size_t i = 0; i < bonds_->keys_.size(); i++) {
            auto obj = json::object();
            obj["key"] = bonds_->keys_[i];
            obj["value"] = bonds_->parameters_[i];
            j["atom"]["data"].push_back(obj);
        }
    }

    std::ofstream f(filename);
    f << j.dump(4) << std::endl;
    f.close();
}

void Parameters::print() const {
    fmt::print("Parameters: {}\n", name_);
    if (common_) {
        fmt::print("Common parameters\n");
        for (size_t i = 0; i < common_->names().size(); i++) {
            fmt::print("{}: {:.3f}\n", common_->names()[i], common_->parameters_[i]);
        }
    }
    if (atoms_) {
        fmt::print("Atom parameters\n");
        for (size_t i = 0; i < atoms_->parameters_.size(); i++) {
            auto &[symbol, cls, type] = atoms_->keys()[i];
            fmt::print("{:2s} {} {}: ", symbol, cls, type);
            for (double val: atoms_->parameters_[i]) {
                fmt::print("{:>-6.3f} ", val);
            }
            fmt::print("\n");
        }
    }
    if (bonds_) {
        fmt::print("Bond parameters\n");
        for (size_t i = 0; i < bonds_->parameters_.size(); i++) {
            auto &[symbol1, symbol2, cls, type] = bonds_->keys()[i];
            fmt::print("{:2s} {:2s} {} {}: ", symbol1, symbol2, cls, type);
            for (double val: bonds_->parameters_[i]) {
                fmt::print("{:>-6.3f} ", val);
            }
            fmt::print("\n");
        }
    }
}


std::function<double(const Atom &)> AtomParameters::parameter(size_t idx) const {

    return [this, idx](const Atom &atom) { return parameters_[atom.atom_type()][idx]; };
}


std::function<double(const Bond &)> BondParameters::parameter(size_t idx) const {

    return [this, idx](const Bond &bond) { return parameters_[bond.bond_type()][idx]; };
}


std::vector<double> Parameters::get_vector() const {
    std::vector<double> parameters;
    if (common_) {
        for (const auto &v: common_->parameters_) {
            parameters.push_back(v);
        }
    }
    if (atoms_) {
        for (const auto &key: atoms_->parameters_) {
            for (const auto &v: key) {
                parameters.push_back(v);
            }
        }
    }
    if (bonds_) {
        for (const auto &key: bonds_->parameters_) {
            for (const auto &v: key) {
                parameters.push_back(v);
            }
        }
    }
    return parameters;
}


void Parameters::set_from_vector(const std::vector<double> &parameters) {
    size_t idx = 0;
    if (common_) {
        for (auto &v: common_->parameters_) {
            v = parameters[idx++];
        }
    }
    if (atoms_) {
        for (auto &key: atoms_->parameters_) {
            for (auto &v: key) {
                v = parameters[idx++];
            }
        }
    }
    if (bonds_) {
        for (auto &key: bonds_->parameters_) {
            for (auto &v: key) {
                v = parameters[idx++];
            }
        }
    }
}


Parameters::Parameters(const MoleculeSet &ms, const std::shared_ptr<Method> &method) {

    method_name_ = method->name();
    name_ = "New parameters";

    auto common_names = method->common_parameters();
    if (!common_names.empty()) {
        std::vector<double> common_values(common_names.size(), 0.0);
        common_ = std::make_unique<CommonParameters>(common_names, common_values);
    }

    auto atom_names = method->atom_parameters();
    if (!atom_names.empty()) {
        auto at = ms.atom_types();
        std::vector<std::vector<double>> atom_values(at.size(), std::vector<double>(atom_names.size()));
        atoms_ = std::make_unique<AtomParameters>(atom_names, atom_values, at);
    }

    auto bond_names = method->bond_parameters();
    if (!bond_names.empty()) {
        auto bt = ms.bond_types();
        std::vector<std::vector<double>> bond_values(bt.size(), std::vector<double>(bond_names.size()));
        bonds_ = std::make_unique<BondParameters>(bond_names, bond_values, bt);
    }
}
