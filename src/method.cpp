//
// Created by krab1k on 6.11.18.
//

#include <vector>
#include <map>
#include <set>
#include <dlfcn.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <Eigen/Dense>
#include <omp.h>

#include "chargefw2.h"
#include "method.h"
#include "parameters.h"
#include "utility/strings.h"


std::vector<RequiredFeatures> Method::get_requirements() const {
    return {};
}


void Method::set_parameters(Parameters *parameters) {
    if (common_parameters_.size() + atom_parameters_.size() + bond_parameters_.size() == 0 and parameters == nullptr) {
        return;
    }
    if (parameters->common() != nullptr and parameters->common()->names() != common_parameters_) {
        fmt::print(stderr, "Parameters don't match\n");
        fmt::print(stderr, "Expected: {}\n", common_parameters_);
        fmt::print(stderr, "Got: {}\n", parameters->common()->names());
        throw std::runtime_error("Invalid common parameters provided");
    }

    if (parameters->atom() != nullptr and parameters->atom()->names() != atom_parameters_) {
        fmt::print(stderr, "Parameters don't match\n");
        fmt::print(stderr, "Expected: {}\n", atom_parameters_);
        fmt::print(stderr, "Got: {}\n", parameters->atom()->names());
        throw std::runtime_error("Invalid atom parameters provided");
    }

    if (parameters->bond() != nullptr and parameters->bond()->names() != bond_parameters_) {
        fmt::print(stderr, "Parameters don't match\n");
        fmt::print(stderr, "Expected: {}\n", bond_parameters_);
        fmt::print(stderr, "Got: {}\n", parameters->bond()->names());
        throw std::runtime_error("Invalid bond parameters provided");
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

bool EEMethod::is_suitable_for_large_molecule() const {
    return true;
}


Eigen::VectorXd EEMethod::solve_EE(const Molecule &molecule,
        const std::function<Eigen::VectorXd(const std::vector<const Atom *> &, double)> &EE_function) const {

    auto method = get_option_value<std::string>("type");
    auto radius = get_option_value<double>("radius");

    if (method != "cover" and molecule.atoms().size() > 80000) {
        fmt::print("Switching to cover as the molecule is too big\n");
        fmt::print("Using radius {}\n", radius);
        method = "cover";
    } else if (method == "full" and molecule.atoms().size() > 20000) {
        fmt::print("Switching to cutoff as the molecule is too big\n");
        fmt::print("Using radius {}\n", radius);
        method = "cutoff";
    }

    if (method == "full") {
        Eigen::setNbThreads(0);
        std::vector<const Atom *> fragment_atoms;
        for (const auto &atom: molecule.atoms()) {
            fragment_atoms.push_back(&atom);
        }

        return EE_function(fragment_atoms, molecule.total_charge());

    } else if (method == "cutoff") {
        const size_t n = molecule.atoms().size();
        Eigen::VectorXd results = Eigen::VectorXd::Zero(n);
        Eigen::setNbThreads(1);

#pragma omp parallel for default(none) shared(results, radius, molecule, EE_function) firstprivate(n)
        for (size_t i = 0; i < n; i++) {
            auto fragment_atoms = molecule.get_close_atoms(molecule.atoms()[i], radius);
            Eigen::VectorXd res = EE_function(fragment_atoms,
                                    static_cast<double>(molecule.total_charge()) * fragment_atoms.size() /
                                    molecule.atoms().size());
            results(i) = res(0);
        }

        double correction = molecule.total_charge() - results.sum();
        correction /= molecule.atoms().size();

        results.array() += correction;
        return results;

    } else /* method == "cover" */ {
        Eigen::setNbThreads(1);

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
        Eigen::VectorXd results = Eigen::VectorXd::Zero(n);
        std::vector<int> charges_count(n, 0);

        std::vector<const Atom *> pivots_vector(pivots.begin(), pivots.end());

#pragma omp parallel for default(none) shared(radius, pivots_vector, neighbors, molecule, charges_count, results, EE_function) firstprivate(n)
        for (size_t i = 0; i < pivots_vector.size(); i++) {
            auto &atom = pivots_vector[i];
            auto fragment_atoms = molecule.get_close_atoms(*atom, radius);
            Eigen::VectorXd res = EE_function(fragment_atoms,
                                    static_cast<double>(molecule.total_charge()) * fragment_atoms.size() / n);

            std::set<size_t> close_atoms = {atom->index()};
            for (const auto &j: neighbors[atom->index()]) {
                close_atoms.insert(j);
                for (const auto &k: neighbors[j]) {
                    close_atoms.insert(k);
                }
            }

            for (const auto &j: close_atoms) {
#pragma omp atomic
                charges_count[j]++;
            }

            for (size_t j = 0; j < fragment_atoms.size(); j++) {
                if (close_atoms.find(fragment_atoms[j]->index()) != close_atoms.end()) {
#pragma omp atomic
                    results(fragment_atoms[j]->index()) += res(j);
                }
            }
        }

        for (long i = 0; i < results.size(); i++) {
            results(i) /= charges_count[i];
        }

        /* 3rd step - correct charges */
        auto correction = (molecule.total_charge() - results.sum()) / n;
        results.array() += correction;
        return results;
    }
}


Method* load_method(const std::string &method_name) {

    std::string file;
    if (ends_with(method_name, ".so")) {
        file = method_name;
    } else {
        file = (std::string(INSTALL_DIR) + "/lib/lib" + method_name + ".so");
    }

    auto handle = dlopen(file.c_str(), RTLD_LAZY);

    auto get_method_handle = (Method *(*)())(dlsym(handle, "get_method"));
    if (!get_method_handle) {
        fmt::print(stderr, "{}\n", dlerror());
        exit(EXIT_FILE_ERROR);
    }

    return (*get_method_handle)();
}
