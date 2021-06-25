//
// Created by krab1k on 05/11/18.
//

#include <fstream>
#include <string>
#include <map>
#include <vector>
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
        throw std::runtime_error("Cannot open file: " + filename);
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
            std::vector<atom_t> keys;
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
            std::vector<bond_t> keys;
            std::vector<std::vector<double>> parameters;
            for (const auto &obj: j["bond"]["data"]) {
                auto key = obj["key"];
                keys.emplace_back(
                        std::make_tuple(key[0].get<std::string>(), key[1].get<std::string>(), key[2].get<std::string>(),
                                        key[3].get<std::string>(), key[4].get<std::string>(), key[5].get<std::string>(),
                                        key[6].get<std::string>(), key[7].get<std::string>()));
                parameters.emplace_back(obj["value"].get<std::vector<double>>());
            }
            bonds_ = std::make_unique<BondParameters>(names, parameters, keys);
        }
    }
    catch (std::exception &) {
        throw std::runtime_error("Incorrect file with parameters: " + filename);
    }
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
            fmt::print("{:2s} {:6s} {:4s}: ", symbol, cls, type);
            for (double val: atoms_->parameters_[i]) {
                fmt::print("{:>-6.3f} ", val);
            }
            fmt::print("\n");
        }
    }
    if (bonds_) {
        fmt::print("Bond parameters\n");
        for (size_t i = 0; i < bonds_->parameters_.size(); i++) {
            auto &[symbol1, cls1, type1, symbol2, cls2, type2, cls_b, type_b] = bonds_->keys()[i];
            fmt::print("{:2s} {:6s} {:4s} {:2s} {:6s} {:4s} {:2s} {:2s}: ", symbol1, cls1, type1, symbol2, cls2, type2, cls_b, type_b);
            for (double val: bonds_->parameters_[i]) {
                fmt::print("{:>-6.3f} ", val);
            }
            fmt::print("\n");
        }
    }
}


std::function<double(const Atom &)> AtomParameters::parameter(size_t idx) const noexcept {

    return [this, idx](const Atom &atom) noexcept { return parameters_[atom.type()][idx]; };
}


std::function<double(const Bond &)> BondParameters::parameter(size_t idx) const noexcept {

    return [this, idx](const Bond &bond) noexcept { return parameters_[bond.type()][idx]; };
}

