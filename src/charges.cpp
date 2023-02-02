//
// Created by krab1k on 8.11.18.
//

#include <string>
#include <fstream>
#include <sstream>
#include <limits>

#include "charges.h"


void Charges::insert(const std::string &name, std::vector<double> charges) {
    names_.push_back(name);
    charges_[name] = std::move(charges);
}


Charges::Charges(const std::string &filename) {
    std::ifstream f(filename);

    if (not f.is_open()) {
        throw std::runtime_error("Cannot open file with charges");
    }

    std::string line;
    while (std::getline(f, line)) {
        std::string name = line;
        names_.push_back(name);

        std::getline(f, line);
        double value = std::numeric_limits<double>::quiet_NaN();
        std::vector<double> values;
        std::stringstream ss(line);
        while (ss >> value) {
            values.push_back(value);
        }
        charges_[name] = values;
    }
    f.close();
}


Charges::Charges(const std::map<std::string, std::vector<double>> &charges)
    : charges_(charges){}
