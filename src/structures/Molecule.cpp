//
// Created by krab1k on 23/10/18.
//

#include "Atom.h"
#include "Molecule.h"
#include <utility>

Molecule::Molecule(QString name, QVector<Atom> atoms, QVector<Bond> bonds) {
    name_ = std::move(name);
    atoms_ = atoms;
    bonds_ = bonds;
}

std::ostream &operator<<(std::ostream &str, const Molecule &molecule) {
    return str << "Molecule: " << molecule.name_.toStdString() << " Atoms: " << molecule.atoms_.size() << " Bonds: "
               << molecule.bonds_.size();
}