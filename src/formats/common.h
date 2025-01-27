#pragma once

#include <string>
#include <set>
#include <gemmi/model.hpp>

#include "../structures/atom.h"

std::string sanitize_name(const std::string &name);

std::string get_unique_name(const std::string &name, const std::set<std::string> &already_used);

std::string get_element_symbol(const std::string &substring);

std::string fix_atom_name(std::string &atom_name);

bool keep_atom(const gemmi::Atom &atom, const gemmi::Residue &residue);
