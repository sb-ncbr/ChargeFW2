//
// Created by krab1k on 29/10/18.
//


#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "Element.h"
#include "PeriodicTable.h"

PeriodicTable::PeriodicTable() {
    std::ifstream file("../share/pte.csv");
    if (!file) {
        std::cerr << "Unable to open periodic table data file data/pte.csv" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    // Read header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream line_stream(line);
        std::string cell;
        std::vector<std::string> cols;

        while(std::getline(line_stream, cell, ',')) {
            cols.emplace_back(cell);
        };

        int index = std::stoi(cols[0]);
        std::string name = cols[1];
        std::string symbol = cols[2];
        float electronegativity = std::stof(cols[11]);
        Element element(index, symbol, name, electronegativity);
        elements_.push_back(element);
        symbol_Z_[symbol] = index;
    }
}

const PeriodicTable &PeriodicTable::pte() {
    static PeriodicTable pte;
    return pte;
}

const Element *PeriodicTable::getElement(const std::string &symbol) const {
    if (!symbol_Z_.count(symbol)) {
        std::cerr << "No such element " << symbol << std::endl;
        exit(EXIT_FAILURE);
    }
    return getElement(symbol_Z_.at(symbol) - 1);
}
