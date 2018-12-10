//
// Created by krab1k on 29/10/18.
//


#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "element.h"
#include "periodic_table.h"
#include "config.h"

PeriodicTable::PeriodicTable() {
    std::string filename(std::string(INSTALL_DIR) + "/share/pte.csv");
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Unable to open periodic table data file: " << filename << std::endl;
        exit(EXIT_INTERNAL_ERROR);
    }

    std::string line;
    // Read header
    std::getline(file, line);

    try {
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
            float covalent_radius = std::stof(cols[10]);
            float vdw_radius = std::stof(cols[11]);
            float electronegativity = std::stof(cols[12]);
            Element element(index, symbol, name, electronegativity, covalent_radius, vdw_radius);
            elements_.push_back(element);
            symbol_Z_[symbol] = index;
        }
    } catch (std::invalid_argument &err) {
        std::cerr << "Unable to read periodic table data file: " << filename << std::endl;
        exit(EXIT_INTERNAL_ERROR);
    }
}

const PeriodicTable &PeriodicTable::pte() {
    static PeriodicTable pte;
    return pte;
}

const Element *PeriodicTable::getElement(const std::string &symbol) const {
    if (!symbol_Z_.count(symbol)) {
        std::cerr << "No such element " << symbol << std::endl;
        exit(EXIT_INTERNAL_ERROR);
    }
    return getElement(symbol_Z_.at(symbol) - 1);
}
