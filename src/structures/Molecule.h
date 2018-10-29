//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <string>
#include <utility>
#include <QString>
#include <QVector>

#include "Atom.h"
#include "Bond.h"

class Molecule {
    QString name_;
    QVector<Atom> atoms_;
    QVector<Bond> bonds_;

public:
    const QVector<Atom> &atoms() const { return atoms_; }

    const QVector<Bond> &bonds() const { return bonds_; }

    const QString &name() const { return name_; }

    Molecule() = default;

    Molecule(QString name, QVector<Atom> atoms, QVector<Bond> bonds);
};
