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
    std::map<std::string, size_t> name_Z_;
public:
    static const PeriodicTable &pte();

    [[nodiscard]] const Element *get_element_by_Z(size_t Z) const { return &elements_[Z]; }

    [[nodiscard]] const Element *get_element_by_symbol(const std::string &symbol) const;

    [[nodiscard]] const Element *get_element_by_name(const std::string &symbol) const;

    PeriodicTable();
};
