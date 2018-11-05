//
// Created by krab1k on 01/11/18.
//

#pragma once


#include <QString>
#include "structures/Atom.h"
#include "structures/MoleculeSet.h"

class Classifier {
public:
    virtual QString get_type(const Atom &atom) const = 0;

    virtual QString name() const = 0;
};

class PlainClassifier : public Classifier {
public:
    QString name() const override { return QString("plain"); }

    QString get_type(const Atom &) const override { return QString("*"); }
};


class HBOClassifier : public Classifier {
public:
    QString name() const override { return QString("hbo"); }

    QString get_type(const Atom &atom) const override;
};