#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <fmt/format.h>
#include <tuple>
#include <vector>
#include <string>

#include "chargefw2.h"
#include "candidates.h"
#include "method.h"
#include "utility/strings.h"


namespace fs = std::filesystem;


std::vector<std::string>
get_valid_parameters(MoleculeSet &ms, bool is_protein, bool permissive_types, const std::string &method_name);


std::vector<fs::path> get_parameter_files() {
    /* Get parameters sorted according to the names and priorities */
    std::vector<fs::path> files;
    for (const auto &set: fs::directory_iterator(fs::path(INSTALL_DIR) / "share/parameters")) {
        files.emplace_back(set.path());
    }
    std::sort(files.begin(), files.end());
    return files;
}


std::vector<std::string>
get_valid_parameters(MoleculeSet &ms, bool is_protein, bool permissive_types, const std::string &method_name) {
    std::vector<std::string> protein_parameters;
    std::vector<std::string> ligand_parameters;

    for (const auto &parameter_file: get_parameter_files()) {
        if (not starts_with(to_lowercase(parameter_file.filename().string()), method_name)) {
            continue;
        }

        auto p = std::make_unique<Parameters>(parameter_file);

        if (method_name != p->method_name()) {
            continue;
        }

        size_t unclassified = ms.classify_set_from_parameters(*p, false, permissive_types);

        auto parameters = parameter_file.filename().string();
        if (!unclassified) {
            if (p->source() == "protein") {
                protein_parameters.emplace_back(parameters);
            } else {
                ligand_parameters.emplace_back(parameters);
            }
        }
    }

    /* Show protein parameters first if the proteins are in the set */
    std::vector<std::string> all_parameters;
    if (is_protein) {
        all_parameters = protein_parameters;
        all_parameters.insert(all_parameters.end(), ligand_parameters.begin(), ligand_parameters.end());
    } else {
        all_parameters = ligand_parameters;
        all_parameters.insert(all_parameters.end(), protein_parameters.begin(), protein_parameters.end());
    }

    return all_parameters;
}


std::vector<std::tuple<std::string, std::vector<std::string>>>
get_suitable_methods(MoleculeSet &ms, bool is_protein, bool permissive_types) {
    std::string filename = (fs::path(INSTALL_DIR) / "share/methods.json").string();
    using json = nlohmann::json;
    json j;
    std::ifstream f(filename);
    if (!f) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    f >> j;
    f.close();

    std::vector<std::tuple<std::string, std::vector<std::string>>> results;

    for (const auto &method_info: j["methods"]) {
        auto method_name = method_info["internal_name"].get<std::string>();
        const auto method = load_method(method_name);

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
            results.emplace_back(std::tuple<std::string, std::vector<std::string>>(method_name, {}));
            continue;
        }

        auto parameters = get_valid_parameters(ms, is_protein, permissive_types, method_name);
        if (not parameters.empty()) {
            results.emplace_back(std::make_tuple(method_name, parameters));
        }
    }

    return results;
}


std::string
best_parameters(MoleculeSet &ms, const Method *method, bool is_protein, bool permissive_types) {
    auto parameters = get_valid_parameters(ms, is_protein, permissive_types, method->internal_name());
    return parameters.empty() ? "" : parameters.front();
}
