#include <string>
#include <vector>

#include "charges.h"


void Charges::insert(const std::string &name, std::vector<double> charges) {
    names_.push_back(name);
    charges_[name] = std::move(charges);
}
