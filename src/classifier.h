//
// Created by krab1k on 01/11/18.
//

#pragma once


#include <string>
#include "structures/atom.h"
#include "structures/molecule_set.h"

class AtomClassifier {
public:
    virtual std::string get_type(const Atom &atom) const = 0;

    virtual std::string name() const = 0;
};


class BondClassifier {
public:
    virtual std::string get_type(const Bond &atom) const = 0;

    virtual std::string name() const = 0;
};


class PlainAtomClassifier : public AtomClassifier {
public:
    std::string name() const override { return std::string("plain"); }

    std::string get_type(const Atom &) const override { return std::string("*"); }
};


class HBOAtomClassifier : public AtomClassifier {
public:
    std::string name() const override { return std::string("hbo"); }

    std::string get_type(const Atom &atom) const override;
};


class PlainBondClassifier : public BondClassifier {
public:
    std::string name() const override { return std::string("plain"); }

    std::string get_type(const Bond &) const override { return std::string("*"); }
};


class BOBondClassifier : public BondClassifier {
public:
    std::string name() const override { return std::string("bo"); }

    std::string get_type(const Bond &bond) const override { return std::to_string(bond.order()); }
};