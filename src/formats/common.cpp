#include <string>
#include <set>
#include <cctype>

#include "common.h"
#include "../config.h"
#include "../structures/atom.h"
#include "../utility/strings.h"


std::string sanitize_name(const std::string &name) {
    std::string new_name;
    for (const auto c: name) {
        if (isalnum(c)) {
            new_name += c;
        } else {
            new_name += '_';
        }
    }
    return new_name;
}


std::string get_unique_name(const std::string &name, const std::set<std::string> &already_used) {
    int count = 0;
    std::string new_name = name;
    while (already_used.contains(new_name)) {
        new_name = name + std::string("_") + std::to_string(count);
        count++;
    }
    return new_name;
}


std::string get_element_symbol(const std::string &substring) {
    auto element_symbol = trim(substring);
    element_symbol = to_lowercase(element_symbol);
    element_symbol[0] = static_cast<char>(std::toupper(element_symbol[0]));

    return element_symbol;
}


std::string fix_atom_name(std::string &atom_name) {

    if (atom_name.front() == '"' and atom_name.back() == '"') {
        return atom_name.substr(1, atom_name.size() - 2);
    }

    return atom_name;
}


bool keep_atom(const gemmi::Atom& atom, const gemmi::Residue& residue) {
    if (not atom.has_altloc() or atom.altloc == 'A') {
        const bool hetatm = residue.het_flag == 'H';
        if (not hetatm or (config::read_hetatm and (residue.name != "HOH" or not config::ignore_water))) {
            return true;
        }
    }
    return false;
}
