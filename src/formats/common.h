//
// Created by krab1k on 31.1.19.
//


#pragma once

#include <string>
#include <vector>

#include "../structures/atom.h"


std::string get_element_symbol(const std::string &substring);

std::string fix_atom_name(std::string &atom_name);

bool is_already_loaded(const std::vector<Atom> &atoms, const std::string &atom_name, int residue_id);
