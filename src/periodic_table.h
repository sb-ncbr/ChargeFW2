//
// Created by krab1k on 29/10/18.
//

#pragma once

#include <string>
#include <vector>
#include <map>

#include "element.h"


class PeriodicTable {
    std::vector<Element> elements_;
    std::map<std::string, size_t> symbol_Z_;
public:
    static const PeriodicTable &pte();

    const Element *getElement(size_t Z) const { return &elements_[Z]; }

    const Element *getElement(const std::string &symbol) const;

    PeriodicTable();
};
