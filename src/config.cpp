#include <fmt/core.h>

#include "chargefw2.h"
#include "config.h"


namespace config {
    std::string mode;
    std::string input_file;
    std::string par_file;
    std::string chg_out_dir;
    std::string log_file;
    std::string method_name;
    bool read_hetatm;
    bool ignore_water;
    bool permissive_types;
}


void check_common_args() {
    if (config::mode != "available-methods" && config::input_file.empty()) {
        fmt::print(stderr, "Input file must be provided\n");
        exit(EXIT_PARAMETER_ERROR);
    }

    if (config::mode == "best-parameters") {
        if (config::method_name.empty()) {
            fmt::print(stderr, "No method selected.\n");
            exit(EXIT_PARAMETER_ERROR);
        }
    }

    if (config::mode == "charges") {
        if (config::chg_out_dir.empty()) {
            fmt::print(stderr, "Directory where to store charges must be provided\n");
            exit(EXIT_PARAMETER_ERROR);
        }
    }
}
