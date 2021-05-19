//
// Created by krab1k on 31.1.19.
//


#pragma once

#include <string>
#include <vector>
#include <set>

#include "../structures/atom.h"

std::string sanitize_name(const std::string &name);

std::string get_unique_name(const std::string &name, const std::set<std::string> &already_used);

std::string get_element_symbol(const std::string &substring);

std::string fix_atom_name(std::string &atom_name);
