//
// Created by krab1k on 31.1.19.
//

#include <string>
#include <boost/algorithm/string.hpp>

#include "common.h"


std::string get_element_symbol(const std::string &substring) {
    auto element_symbol = boost::trim_copy(substring);
    boost::to_lower(element_symbol);
    element_symbol[0] = static_cast<char>(std::toupper(element_symbol[0]));

    return element_symbol;
}


std::string fix_atom_name(std::string &atom_name) {

    if (atom_name.front() == '"' and atom_name.back() == '"') {
        return atom_name.substr(1, atom_name.size() - 2);
    }

    return atom_name;
}
