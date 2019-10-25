//
// Created by krab1k on 29/10/18.
//


#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <fmt/format.h>

#include "chargefw2.h"
#include "element.h"
#include "periodic_table.h"


PeriodicTable::PeriodicTable() {
    std::string filename(std::string(INSTALL_DIR) + "/share/pte.csv");
    std::ifstream file(filename);
    if (!file) {
        fmt::print(stderr, "Unable to open periodic table data file: {}\n", filename);
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

            while (std::getline(line_stream, cell, ',')) {
                cols.emplace_back(cell);
            }

            size_t index = std::stoul(cols[0]);
            std::string name = cols[1];
            std::string symbol = cols[2];

            int period = std::stoi(cols[4]);
            int group = std::stoi(cols[5]);
            // Change units from pm to A
            double covalent_radius = std::stod(cols[10]) / 100.0;
            double vdw_radius = std::stod(cols[11]) / 100.0;
            double electronegativity = std::stod(cols[12]);
            double electron_affinity = std::stod(cols[14]);
            double ionization_potential = std::stod(cols[13]);
            Element element(index, symbol, name, electronegativity, covalent_radius, vdw_radius, period, group,
                            electron_affinity, ionization_potential);
            elements_.push_back(element);
            symbol_Z_[symbol] = index;
            name_Z_[name] = index;
        }
    } catch (std::invalid_argument &) {
        fmt::print(stderr, "Unable to read periodic table data file: {}\n", filename);
        exit(EXIT_INTERNAL_ERROR);
    }
}


const PeriodicTable &PeriodicTable::pte() {
    static PeriodicTable pte;
    return pte;
}


const Element *PeriodicTable::get_element_by_symbol(const std::string &symbol) const {
    /* Treat deuterium as hydrogen */
    if (symbol == "D") {
        return get_element_by_Z(0);
    }

    if (!symbol_Z_.count(symbol)) {
        throw std::runtime_error(fmt::format("No such element: {}", symbol));
    }
    return get_element_by_Z(symbol_Z_.at(symbol) - 1);
}


const Element *PeriodicTable::get_element_by_name(const std::string &name) const {
    if (!name_Z_.count(name)) {
        throw std::runtime_error(fmt::format("No such element: {}", name));
    }
    return get_element_by_Z(name_Z_.at(name) - 1);
}
