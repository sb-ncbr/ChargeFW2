#pragma once

#include <boost/program_options.hpp>

#include "method.h"


void setup_method_options(Method *method, const boost::program_options::parsed_options& parsed);

boost::program_options::parsed_options parse_args(int argc, char **argv);
