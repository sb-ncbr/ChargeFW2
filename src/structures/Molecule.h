//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <iostream>
#include <utility>
#include <QString>
#include <QVector>
#include <vector>
#include <tuple>
#include <memory>

#include "Atom.h"
#include "Bond.h"

class Molecule {
    QString name_;
    std::unique_ptr<std::vector<Atom> > atoms_;
    std::unique_ptr<std::vector<Bond> > bonds_;

public:
    const std::vector<Atom> &atoms() const { return *atoms_; }

    const std::vector<Bond> &bonds() const { return *bonds_; }

    const QString &name() const { return name_; }

    Molecule() = default;

    Molecule(QString name, std::unique_ptr<std::vector<Atom> > atoms, std::unique_ptr<std::vector<Bond> > bonds);

    friend std::ostream &operator<<(std::ostream &str, const Molecule &molecule);

    friend class MoleculeSet;
};
