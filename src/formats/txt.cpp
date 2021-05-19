//
// Created by krab1k on 30.1.19.
//

#include <fstream>
#include <string>
#include <fmt/format.h>

#include "chargefw2.h"
#include "../charges.h"
#include "txt.h"
#include "../utility/utility.h"
#include "../utility/strings.h"


void TXT::save_charges(const MoleculeSet &, const Charges &charges, const std::string &filename) {
    std::ofstream file(filename);
    if (!file) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    for (const auto &name: charges.names()) {
        file << to_uppercase(name) << std::endl;
        file << charges[name] << std::endl;
    }

    file.close();
}
