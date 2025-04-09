#include <boost/program_options.hpp>
#include <fmt/core.h>
#include <fmt/ostream.h>

#include "chargefw2.h"
#include "config.h"
#include "options.h"


boost::program_options::parsed_options parse_args(int argc, char **argv) {
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "Prints this help")
            ("mode", po::value<std::string>()->required(), "Mode")
            ("input-file", po::value<std::string>()->default_value(""), "Input file")
            ("par-file", po::value<std::string>()->default_value(""), "File with parameters (json)")
            ("chg-out-dir", po::value<std::string>()->default_value(""), "Directory to output charges to")
            ("log-file", po::value<std::string>()->default_value(""), "Log file")
            ("read-hetatm", po::bool_switch()->default_value(false), "Read HETATM records from PDB/mmCIF files")
            ("ignore-water", po::bool_switch()->default_value(false), "Discard water molecules from PDB/mmCIF files")
            ("permissive-types", po::bool_switch()->default_value(false), "Use similar parameters for similar atom/bond types if no exact match is found")
            ("method", po::value<std::string>()->default_value(""), "Method");

    try {
        po::variables_map vm;
        po::parsed_options parsed = po::command_line_parser(argc, argv)
                .options(desc)
                .style(po::command_line_style::unix_style ^ po::command_line_style::allow_guessing)
                .allow_unregistered()
                .run();

        po::store(parsed, vm);
        if (vm.contains("help")) {
            fmt::print("ChargeFW2 (version {})\n", VERSION);
            fmt::print("by Tomáš Raček (2018-2023)\n");

#if FMT_VERSION > 90000
            fmt::print("{}", fmt::streamed(desc));
#else
            fmt::print("{}", desc);
#endif
            exit(EXIT_SUCCESS);
        }
        po::notify(vm);

        config::mode = vm["mode"].as<std::string>();
        config::input_file = vm["input-file"].as<std::string>();
        config::par_file = vm["par-file"].as<std::string>();
        config::chg_out_dir = vm["chg-out-dir"].as<std::string>();
        config::log_file = vm["log-file"].as<std::string>();
        config::method_name = vm["method"].as<std::string>();
        config::read_hetatm = vm["read-hetatm"].as<bool>();
        config::ignore_water = vm["ignore-water"].as<bool>();
        config::permissive_types = vm["permissive-types"].as<bool>();

        return parsed;
    } catch (const std::exception &e) {
        fmt::print(stderr, "Incorrect arguments: {}\n", e.what());
        exit(EXIT_PARAMETER_ERROR);
    }
}


void setup_method_options(Method *method, const boost::program_options::parsed_options &parsed) {
    namespace po = boost::program_options;

    po::options_description method_options("Method options");

    for (const auto &[opt, info]: method->get_options()) {
        if (!info.choices.empty()) {
            if (std::ranges::find(info.choices, info.default_value) == info.choices.end()) {
                fmt::print(stderr, "Default value: {} not in possible choices\n", info.default_value);
                exit(EXIT_INTERNAL_ERROR);
            }
        }
        method_options.add_options()(std::string("method-" + opt).c_str(), po::value<std::string>(),
                                     info.description.c_str());
    }

    po::variables_map vm;
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    try {
        po::store(po::command_line_parser(opts).options(method_options).run(), vm);
    } catch (std::exception &e) {
        fmt::print(stderr, "Incorrect arguments: {}\n", e.what());
        exit(EXIT_PARAMETER_ERROR);
    }

    for (const auto &[opt, info]: method->get_options()) {
        std::string opt_name = std::string("method-" + opt);
        if (vm.contains(opt_name)) {
            std::string val = vm[opt_name].as<std::string>();
            if (!info.choices.empty()) {
                if (std::ranges::find(info.choices, val) == info.choices.end()) {
                    fmt::print(stderr, "Provided value: {} not in possible choices\n", val);
                    exit(EXIT_INTERNAL_ERROR);
                }
            }
            method->set_option_value(opt, val);
        } else {
            method->set_option_value(opt, info.default_value);
        }
    }
}
