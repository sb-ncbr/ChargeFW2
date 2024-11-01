#include <string>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fmt/os.h>

#include "chargefw2.h"
#include "../charges.h"
#include "txt.h"
#include "../utility/strings.h"


void TXT::save_charges(const MoleculeSet &, const Charges &charges, const std::string &filename) {
    try {
        auto file = fmt::output_file(filename);
        for (const auto &name: charges.names()) {
            file.print("{}\n", to_uppercase(name));
            file.print("{:.5f}\n", fmt::join(charges[name], " "));
        }
        file.close();
    } catch (std::system_error &e) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }
}
