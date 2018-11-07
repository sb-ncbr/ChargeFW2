//
// Created by krab1k on 23/10/18.
//

#include "Atom.h"
#include "Bond.h"
#include "Molecule.h"
#include <utility>
#include <tuple>
#include <string>
#include <vector>
#include <iostream>


Molecule::Molecule(std::string name, std::unique_ptr<std::vector<Atom> > atoms, std::unique_ptr<std::vector<Bond> > bonds) {
    name_ = std::move(name);
    atoms_ = std::move(atoms);
    bonds_ = std::move(bonds);
}

std::ostream &operator<<(std::ostream &str, const Molecule &molecule) {
    return str << "Molecule: " << molecule.name_ << " Atoms: " << molecule.atoms_->size() << " Bonds: "
               << molecule.bonds_->size();
}
