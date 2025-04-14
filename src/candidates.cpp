#include <fmt/core.h>
#include <iterator>
#include <memory>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fmt/format.h>
#include <optional>
#include <tuple>
#include <vector>
#include <string>

#include "chargefw2.h"
#include "candidates.h"
#include "exceptions/parameter_exception.h"
#include "method.h"
#include "parameters.h"
#include "utility/strings.h"
#include "utility/install.h"


namespace fs = std::filesystem;


std::vector<std::unique_ptr<Parameters>>
get_valid_parameters(MoleculeSet &ms, bool is_protein, bool permissive_types, const std::string &method_name);


std::vector<fs::path> get_parameter_files() {
    /* Get parameters sorted according to the names and priorities */
    std::vector<fs::path> files;
    for (const auto &set: fs::directory_iterator(InstallPaths::datadir() / "parameters")) {
        files.emplace_back(set.path());
    }
    std::ranges::sort(files);
    return files;
}


std::vector<std::unique_ptr<Parameters>>
get_valid_parameters(MoleculeSet &ms, bool is_protein, bool permissive_types, const std::string &method_name) {
    std::vector<std::unique_ptr<Parameters>> protein_parameters;
    std::vector<std::unique_ptr<Parameters>> ligand_parameters;

    for (const auto &parameter_file: get_parameter_files()) {
        if (not to_lowercase(parameter_file.filename().string()).starts_with(method_name)) {
            continue;
        }

        auto p = std::make_unique<Parameters>(parameter_file);

        if (method_name != p->method_name()) {
            continue;
        }

        size_t unclassified = ms.classify_set_from_parameters(*p, false, permissive_types);

        if (!unclassified) {
            if (p->source() == "protein") {
                protein_parameters.emplace_back(std::move(p));
            } else {
                ligand_parameters.emplace_back(std::move(p));
            }
        }
    }

    /* Show protein parameters first if the proteins are in the set */
    std::vector<std::unique_ptr<Parameters>> all_parameters;
    if (is_protein) {
        all_parameters = std::move(protein_parameters);
        all_parameters.insert(
            all_parameters.end(), 
            std::make_move_iterator(ligand_parameters.begin()), 
            std::make_move_iterator(ligand_parameters.end()));
    } else {
        all_parameters = std::move(ligand_parameters);
        all_parameters.insert(
            all_parameters.end(), 
            std::make_move_iterator(protein_parameters.begin()), 
            std::make_move_iterator(protein_parameters.end()));
    }

    return all_parameters;
}


std::vector<std::tuple<Method*, std::vector<std::unique_ptr<Parameters>>>>
get_suitable_methods(MoleculeSet &ms, bool is_protein, bool permissive_types) {
    std::vector<std::tuple<Method*, std::vector<std::unique_ptr<Parameters>>>> results;

    auto methods = get_available_methods();
    std::sort(methods.begin(), methods.end(), [](const auto &a, const auto &b) {
        return a->metadata().priority > b->metadata().priority;
    });

    for (auto &method: methods) {
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
            results.emplace_back(method, std::vector<std::unique_ptr<Parameters>>{});
            continue;
        }

        auto parameters = get_valid_parameters(ms, is_protein, permissive_types, method->metadata().internal_name);
        if (not parameters.empty()) {
            results.emplace_back(method, std::move(parameters));
        }
    }

    return results;
}


std::optional<std::unique_ptr<Parameters>>
best_parameters(MoleculeSet &ms, const Method *method, bool is_protein, bool permissive_types) {
    if (not method->has_parameters()) {
        throw ParameterException("Method uses no parameters");
    }
    
    auto parameters = get_valid_parameters(ms, is_protein, permissive_types, method->metadata().internal_name);
    
    if (parameters.empty()) {
        return std::nullopt;
    }

    return std::move(parameters.front());
}
