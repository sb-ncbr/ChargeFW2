#pragma once

#include <memory>
#include <tuple>

#include "structures/molecule_set.h"
#include "method.h"


std::vector<std::filesystem::path> get_parameter_files();

std::vector<std::tuple<const Method *, std::vector<std::unique_ptr<Parameters>>>>
get_suitable_methods(MoleculeSet& ms, const std::vector<std::unique_ptr<Method>>& methods, bool is_protein, bool permissive_types = false);

std::optional<std::unique_ptr<Parameters>>
best_parameters(MoleculeSet& ms, const Method& method, bool is_protein, bool permissive_types = false);
