#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <boost/algorithm/string.hpp>

#include "chargefw2.h"
#include "candidates.h"
#include "method.h"
#include "utility/utility.h"


std::vector<std::string> get_parameter_files();


std::vector<std::string> get_parameter_files() {
    /* Get parameters sorted according to the names and priorities */
    std::vector<std::string> files;
    for (const auto &set: std::filesystem::directory_iterator(std::string(INSTALL_DIR) + "/share/parameters")) {
        files.emplace_back(set.path().string());
    }
    std::sort(files.begin(), files.end());
    return files;
}


void get_suitable_methods(MoleculeSet &ms, bool is_protein, bool permissive_types) {
    std::string filename = std::string(INSTALL_DIR) + "/share/methods.json";
    using json = nlohmann::json;
    json j;
    std::ifstream f(filename);
    if (!f) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    f >> j;
    f.close();

    for (const auto &method_info: j["methods"]) {
        auto method_name = method_info["internal_name"].get<std::string>();
        auto method = load_method(method_name);

        bool suitable = true;
        for (const auto &molecule: ms.molecules()) {
            if (not method->is_suitable_for_molecule(molecule) or
                (molecule.atoms().size() > LARGE_MOLECULE_ATOM_COUNT and
                 not method->is_suitable_for_large_molecule())) {
                suitable = false;
                break;
            }
        }

        if (not suitable) {
            continue;
        }

        /* Methods without parameters should be suitable */
        if (not method->has_parameters()) {
            fmt::print("{}\n", method_name);
            continue;
        }

        bool parameters_found = false;
        for (const auto &set: get_parameter_files()) {
            if (not boost::starts_with(to_lowercase(std::filesystem::path(set).filename().string()), method_name)) {
                continue;
            }

            auto p = std::make_unique<Parameters>(set);

            if (method_name != p->method_name()) {
                continue;
            }

            if ((is_protein and p->source() != "protein") or (not is_protein and p->source() == "protein")) {
                continue;
            }

            size_t unclassified = ms.classify_set_from_parameters(*p, false, permissive_types);

            if (!unclassified) {
                if (not parameters_found) {
                    fmt::print("{} ", method_name);
                    parameters_found = true;
                }
                fmt::print("{} ", std::filesystem::path(set).filename().string());
            }
        }
        if (parameters_found) {
            fmt::print("\n");
        }
    }
}


std::string
best_parameters(MoleculeSet &ms, const std::shared_ptr<Method> &method, bool is_protein, bool permissive_types) {
    std::string best_name;
    std::string best_name_permissive;
    size_t best_unclassified = ms.molecules().size();
    size_t best_unclassified_permissive = ms.molecules().size();

    auto internal = method->internal_name();
    for (const auto &set: get_parameter_files()) {
        if (not boost::starts_with(to_lowercase(std::filesystem::path(set).filename().string()), internal)) {
            continue;
        }

        auto p = std::make_unique<Parameters>(set);

        if (internal != p->method_name()) {
            continue;
        }

        if ((is_protein and p->source() != "protein") or (not is_protein and p->source() == "protein")) {
            continue;
        }

        size_t unclassified = ms.classify_set_from_parameters(*p, false);
        size_t unclassified_permissive = ms.molecules().size();

        if (unclassified != 0 and permissive_types) {
            unclassified_permissive = ms.classify_set_from_parameters(*p, false, permissive_types);
        }

        // If all molecules are covered by the parameters, we found our best
        if (!unclassified) {
            return set;
        }

        if (unclassified < best_unclassified) {
            best_unclassified = unclassified;
            best_name = set;
        }

        if (permissive_types and unclassified_permissive < best_unclassified_permissive) {
            best_unclassified_permissive = unclassified_permissive;
            best_name_permissive = set;
        }
    }

    if (permissive_types and best_unclassified_permissive < best_unclassified) {
        return best_name_permissive;
    } else {
        return best_name;
    }
}
