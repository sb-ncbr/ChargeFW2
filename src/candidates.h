#pragma once

#include "structures/molecule_set.h"
#include "method.h"

std::vector<std::filesystem::path> get_parameter_files();

std::vector<std::tuple<MethodMetadata, std::vector<ParametersMetadata>>>
get_suitable_methods(MoleculeSet &ms, bool is_protein, bool permissive_types = false);

std::optional<ParametersMetadata>
best_parameters(MoleculeSet &ms, const Method *method, bool is_protein, bool permissive_types = false);
