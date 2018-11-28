//
// Created by krab1k on 01/11/18.
//

#pragma once


#include <string>
#include "structures/atom.h"
#include "structures/molecule_set.h"

class Classifier {
public:
    virtual std::string get_type(const Atom &atom) const = 0;

    virtual std::string name() const = 0;
};

class PlainClassifier : public Classifier {
public:
    std::string name() const override { return std::string("plain"); }

    std::string get_type(const Atom &) const override { return std::string("*"); }
};


class HBOClassifier : public Classifier {
public:
    std::string name() const override { return std::string("hbo"); }

    std::string get_type(const Atom &atom) const override;
};