//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "Atom.h"
#include "Bond.h"

class Molecule {
    std::string name_;
    std::vector<Atom> atoms_;
    std::vector<Bond> bonds_;

public:
    Molecule(std::string name, std::vector<Atom> atoms, std::vector<Bond> bonds) : name_{std::move(name)},
                                                                                   atoms_{std::move(atoms)},
                                                                                   bonds_{std::move(bonds)} {};

    std::vector<Atom> &atoms() {return atoms_;}
};
