#include <string>
#include <format>
#include <fstream>

#include "../charges.h"
#include "txt.h"
#include "../utility/strings.h"
#include "../utility/exceptions.h"


void TXT::save_charges(const MoleculeSet &, const Charges &charges, const std::string &filename) {
    std::ofstream out(filename);
    if (!out) {
        throw FileException(std::format("Cannot open file: {}", filename));
    }

    for (const auto &name: charges.names()) {
        std::println(out, "{}", to_uppercase(name));
        auto const& vals = charges[name];

        std::string_view sep{};
        for (double v : vals) {
            std::print(out, "{}{:.5f}", sep, v);
            sep = " ";
        }
        std::print(out, "\n");
    }
}
