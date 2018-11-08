//
// Created by krab1k on 8.11.18.
//

#include <fstream>

#include "Charges.h"
#include "utility/Utility.h"

void Charges::insert(const std::string &name, std::vector<double> charges) {
    names_.push_back(name);
    charges_[name] = std::move(charges);
}

void Charges::save_to_file(const std::string &filename) {

    std::ofstream f(filename);
    for(const auto &name: names_) {
        f << name << std::endl;
        f << charges_[name] << std::endl;
    }

    f.close();
}
