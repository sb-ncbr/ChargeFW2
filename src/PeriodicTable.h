//
// Created by krab1k on 29/10/18.
//

#pragma once

#include <string>
#include <vector>
#include <map>

#include "Element.h"

class PeriodicTable {
    std::vector<Element> elements_;
    std::map<std::string, int> symbol_Z_;
public:
    static const PeriodicTable &pte();

    const Element *getElement(int Z) const { return &elements_[Z]; }

    const Element *getElement(const std::string &symbol) const;

    PeriodicTable();
};