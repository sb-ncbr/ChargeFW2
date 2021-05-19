//
// Created by krab1k on 28/11/18.
//

#pragma once

#include <string>
#include <boost/program_options.hpp>

#include "method.h"


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


boost::program_options::parsed_options parse_args(int argc, char **argv);

void check_common_args();

void setup_method_options(Method *method, const boost::program_options::parsed_options& parsed);
