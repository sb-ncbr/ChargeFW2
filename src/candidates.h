#pragma once

#include "structures/molecule_set.h"
#include "method.h"
#include <memory>

std::vector<std::filesystem::path> get_parameter_files();

std::vector<std::tuple<Method*, std::vector<std::unique_ptr<Parameters>>>>
get_suitable_methods(MoleculeSet &ms, bool is_protein, bool permissive_types = false);

std::optional<std::unique_ptr<Parameters>>
best_parameters(MoleculeSet &ms, const Method *method, bool is_protein, bool permissive_types = false);
