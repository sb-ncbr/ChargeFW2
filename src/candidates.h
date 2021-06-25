//
// Created by krab1k on 8.4.19.
//

#pragma once

#include "structures/molecule_set.h"
#include "method.h"

std::vector<std::filesystem::path> get_parameter_files();

std::vector<std::tuple<std::string, std::vector<std::string>>>
get_suitable_methods(MoleculeSet &ms, bool is_protein, bool permissive_types = false);

std::string
best_parameters(MoleculeSet &ms, const Method *method, bool is_protein, bool permissive_types = false);
