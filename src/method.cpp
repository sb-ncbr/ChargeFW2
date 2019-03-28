//
// Created by krab1k on 6.11.18.
//

#include <vector>
#include <map>
#include <set>
#include <numeric>
#include <fmt/format.h>
#include <fmt/ranges.h>

#include "chargefw2.h"
#include "method.h"
#include "parameters.h"
#include "utility/utility.h"


std::vector<RequiredFeatures> Method::get_requirements() const {
    return {};
}


void Method::set_parameters(Parameters *parameters) {
    if (common_parameters_.size() + atom_parameters_.size() + bond_parameters_.size() == 0 and parameters == nullptr) {
        return;
    }
    if (parameters->common() != nullptr and parameters->common()->names() != common_parameters_) {
        fmt::print(stderr, "Invalid common parameters provided\n");
        fmt::print(stderr, "Expected: {}\n", common_parameters_);
        fmt::print(stderr, "Got: {}\n", parameters->common()->names());
        exit(EXIT_FILE_ERROR);
    }

    if (parameters->atom() != nullptr and parameters->atom()->names() != atom_parameters_) {
        fmt::print(stderr, "Invalid atom parameters provided\n");
        fmt::print(stderr, "Expected: {}\n", atom_parameters_);
        fmt::print(stderr, "Got: {}\n", parameters->atom()->names());
        exit(EXIT_FILE_ERROR);
    }

    if (parameters->bond() != nullptr and parameters->bond()->names() != bond_parameters_) {
        fmt::print(stderr, "Invalid bond parameters provided\n");
        fmt::print(stderr, "Expected: {}\n", bond_parameters_);
        fmt::print(stderr, "Got: {}\n", parameters->bond()->names());
        exit(EXIT_FILE_ERROR);
    }
    parameters_ = parameters;
}


template<>
std::string Method::get_option_value<std::string>(const std::string &name) const {
    return option_values_.at(name);
}


bool Method::is_suitable_for_molecule(const Molecule &) const {
    return true;
}


std::string Method::internal_name() const {
    auto name = to_lowercase(name_);
    name.erase(std::remove_if(name.begin(), name.end(), [](char c) { return !std::isalnum(c); }), name.end());
    return name;
}


template<>
double Method::get_option_value<double>(const std::string &name) const {
    return std::stod(option_values_.at(name));
}


template<>
int Method::get_option_value<int>(const std::string &name) const {
    return std::stoi(option_values_.at(name));
}


std::vector<double> EEMethod::calculate_charges(const Molecule &molecule) const {

    auto method = get_option_value<std::string>("type");
    auto radius = get_option_value<double>("radius");
    std::vector<const Atom *> fragment_atoms;

    if (method == "full" and molecule.atoms().size() > 50000) {
        fmt::print(stderr, "Switching to cutoff as the molecule is too big\n");
        fmt::print(stderr, "Using radius {}\n", radius);
        method = "cover";
    } else if (method == "full" and molecule.atoms().size() > 20000) {
        fmt::print(stderr, "Switching to cover as the molecule is too big\n");
        fmt::print(stderr, "Using radius {}\n", radius);
        method = "cutoff";
    }

    if (method == "full") {
        for (const auto &atom: molecule.atoms()) {
            fragment_atoms.push_back(&atom);
        }

        return solve_system(fragment_atoms, molecule.total_charge());

    } else if (method == "cutoff") {
        std::vector<double> results;
        for (const auto &atom: molecule.atoms()) {
            fragment_atoms = molecule.get_close_atoms(atom, radius);
            auto res = solve_system(fragment_atoms,
                                    static_cast<double>(molecule.total_charge()) * fragment_atoms.size() /
                                    molecule.atoms().size());
            results.push_back(res[0]);
        }

        double correction = molecule.total_charge();
        for (auto val: results) {
            correction -= val;
        }

        correction /= molecule.atoms().size();

        for (auto &val: results) {
            val += correction;
        }

        return results;
    } else /* method == "cover" */ {
        const size_t n = molecule.atoms().size();

        /* 1st step - identify pivots */
        std::map<size_t, std::set<size_t>> neighbors;
        for (const auto &bond: molecule.bonds()) {
            neighbors[bond.first().index()].insert(bond.second().index());
            neighbors[bond.second().index()].insert(bond.first().index());
        }

        std::set<size_t> all;
        for (size_t i = 0; i < n; i++) {
            all.insert(i);
        }

        std::map<size_t, std::set<size_t>> bonding_sizes;
        for (const auto &[key, val]: neighbors) {
            bonding_sizes[val.size()].insert(key);
        }

        std::set<const Atom *> pivots;
        for (auto it = bonding_sizes.rbegin(); it != bonding_sizes.rend(); it++) {
            for (const auto &idx: it->second) {
                if (all.find(idx) != all.end()) {
                    pivots.insert(&molecule.atoms()[idx]);
                    for (const auto &neighbor: neighbors[idx]) {
                        all.erase(neighbor);
                    }
                }
            }
        }

        /* 2nd step - solve EEM for fragments, sum up charges */
        std::vector<double> results(n, 0);
        std::vector<int> charges_count(n, 0);

        for (const auto &atom: pivots) {
            fragment_atoms = molecule.get_close_atoms(*atom, radius);
            auto res = solve_system(fragment_atoms,
                                    static_cast<double>(molecule.total_charge()) * fragment_atoms.size() / n);

            std::set<size_t> close_atoms = {atom->index()};
            for (const auto &i: neighbors[atom->index()]) {
                close_atoms.insert(i);
                for (const auto &j: neighbors[i]) {
                    close_atoms.insert(j);
                }
            }

            for (const auto &i: close_atoms) {
                charges_count[i]++;
            }

            for (size_t i = 0; i < fragment_atoms.size(); i++) {
                if (close_atoms.find(fragment_atoms[i]->index()) != close_atoms.end()) {
                    results[fragment_atoms[i]->index()] += res[i];
                }
            }
        }

        for (size_t i = 0; i < results.size(); i++) {
            results[i] /= charges_count[i];
        }

        /* 3rd step - correct charges */
        auto sum = std::accumulate(results.begin(), results.end(), 0.0);
        auto correction = (molecule.total_charge() - sum) / n;

        for (auto &r: results) {
            r += correction;
        }

        return results;
    }
}
