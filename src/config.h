//
// Created by krab1k on 28/11/18.
//

#pragma once

#include <string>

namespace config {
    extern std::string mode;
    extern std::string input_file;
    extern std::string par_file;
    extern std::string ref_chg_file;
    extern std::string chg_out_dir;
    extern std::string method_name;
    extern std::string log_file;
    extern bool read_hetatm;
    extern bool ignore_water;
    extern bool permissive_types;
}


void check_common_args();
