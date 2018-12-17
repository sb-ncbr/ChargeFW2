#include <iostream>
#include <boost/program_options.hpp>
#include <boost/dll/import.hpp>
#include <memory>
#include <ctime>
#include <algorithm>
#include <iomanip>
#include <filesystem>

#include "formats/sdf.h"
#include "structures/molecule_set.h"
#include "parameters.h"
#include "classifier.h"
#include "charges.h"
#include "method.h"
#include "config.h"

std::string best_parameters(const MoleculeSet &ms, const boost::shared_ptr<Method> &method);

int main(int argc, char **argv) {
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "Prints this help")
            ("mode", po::value<std::string>()->required(), "Mode")
            ("sdf-file", po::value<std::string>(), "Input SDF file")
            ("par-file", po::value<std::string>(), "File with parameters (json)")
            ("chg-file", po::value<std::string>(), "File to output charges to")
            ("method", po::value<std::string>(), "Method");

    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(argc, argv)
            .options(desc)
            .style(po::command_line_style::unix_style ^ po::command_line_style::allow_guessing)
            .allow_unregistered()
            .run();
    try {

        po::store(parsed, vm);

        if (vm.count("help")) {
            std::cout << "ChargeFW2 (version " << VERSION << ")" << std::endl;
            std::cout << "by Tomáš Raček (2018)\n" << std::endl;
            std::cout << desc << std::endl;
            exit(EXIT_SUCCESS);
        }
        po::notify(vm);
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        exit(EXIT_PARAMETER_ERROR);
    }

    if (!vm.count("sdf-file")) {
        std::cerr << "SDF must be provided" << std::endl;
        exit(EXIT_PARAMETER_ERROR);
    }

    auto sdf_name = vm["sdf-file"].as<std::string>();
    SDF reader;
    MoleculeSet m = reader.read_file(sdf_name);

    auto mode = vm["mode"].as<std::string>();
    if (mode == "info") {
        auto hbo = HBOAtomClassifier();
        m.classify_atoms(hbo);
        m.info();

    } else if (mode == "charges") {
        if (!vm.count("chg-file")) {
            std::cerr << "File where to store charges must be provided" << std::endl;
            exit(EXIT_PARAMETER_ERROR);
        }

        if (!vm.count("method")) {
            std::cerr << "No method selected" << std::endl;
            exit(EXIT_PARAMETER_ERROR);
        }

        auto chg_name = vm["chg-file"].as<std::string>();
        auto method_name = vm["method"].as<std::string>();

        boost::shared_ptr<Method> method;

        try {
            method = boost::dll::import<Method>(std::string(INSTALL_DIR) + "/lib/" + method_name, "method",
                                                boost::dll::load_mode::append_decorations);
        } catch (std::exception &e) {
            std::cerr << "Unable to load method " << method_name << std::endl;
            exit(EXIT_PARAMETER_ERROR);
        }

        po::options_description method_options("Method options");
        for (const auto &[opt, info]: method->get_options()) {
            if(!info.choices.empty()) {
                if (std::find(info.choices.begin(), info.choices.end(), info.default_value) == info.choices.end()) {
                    std::cerr << "Default value: " << info.default_value << " not in possible choices" << std::endl;
                    exit(EXIT_INTERNAL_ERROR);
                }
            }
            method_options.add_options()(std::string("method-" + opt).c_str(), po::value<std::string>(),
                                         info.description.c_str());
        }

        std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        po::store(po::command_line_parser(opts).options(method_options).run(), vm);

        for (const auto &[opt, info]: method->get_options()) {
            std::string opt_name = std::string("method-" + opt);
            if (vm.count(opt_name)){
                std::string val = vm[opt_name].as<std::string>();
                if(!info.choices.empty()) {
                    if (std::find(info.choices.begin(), info.choices.end(), val) == info.choices.end()) {
                        std::cerr << "Provided value: " << val << " not in possible choices" << std::endl;
                        exit(EXIT_INTERNAL_ERROR);
                    }
                }
                method->set_option_value(opt, val);
            }
            else {
                method->set_option_value(opt, info.default_value);
            }
        }

        auto p = std::unique_ptr<Parameters>();

        m.classify_atoms(PlainAtomClassifier());
        m.info();
        std::cout << std::endl;

        if (method->has_parameters()) {
            std::string par_name;
            if (!vm.count("par-file")) {
                par_name = best_parameters(m, method);
                std::cout << "Best parameters found: " << par_name << std::endl;
            } else {
                par_name = vm["par-file"].as<std::string>();
            }

            p = std::make_unique<Parameters>(par_name);

            std::cout << "Parameters:" << std::endl;
            p->print();

            m.classify_set_from_parameters(*p);

            method->set_parameters(p.get());
        } else {
            auto plain = PlainAtomClassifier();
            m.classify_atoms(plain);
        }

        auto charges = Charges();

        clock_t begin = clock();

        for (auto &mol: m.molecules()) {
            charges.insert(mol.name(), method->calculate_charges(mol));
        }

        clock_t end = clock();

        std::cout << "Computation took " << std::setprecision(2) << double(end - begin) / CLOCKS_PER_SEC << " seconds"
                  << std::endl;

        charges.save_to_file(chg_name);

    } else if (mode == "best-parameters") {
        if (!vm.count("method")) {
            std::cerr << "No method selected" << std::endl;
            exit(EXIT_PARAMETER_ERROR);
        }

        auto method_name = vm["method"].as<std::string>();

        boost::shared_ptr<Method> method;

        try {
            method = boost::dll::import<Method>(std::string(INSTALL_DIR) + "/lib/" + method_name, "method",
                                                boost::dll::load_mode::append_decorations);
        } catch (std::exception &e) {
            std::cerr << "Unable to load method " << method_name << std::endl;
            exit(EXIT_PARAMETER_ERROR);
        }

        if (!method->has_parameters()) {
            std::cerr << "Method uses no parameters" << std::endl;
            exit(EXIT_PARAMETER_ERROR);
        }

        std::cout << "Best parameters are: " << best_parameters(m, method) << std::endl;

    } else {
        std::cerr << "Unknown mode: " << mode << std::endl;
        exit(EXIT_PARAMETER_ERROR);
    }

    return 0;
}

std::string best_parameters(const MoleculeSet &ms, const boost::shared_ptr<Method> &method) {
    std::map<std::string, int> missing;

    for (const auto &set: std::filesystem::directory_iterator(std::string(INSTALL_DIR) + "/share/parameters")) {
        auto p = std::make_unique<Parameters>(set.path());

        if (method->name() != p->method_name())
            continue;

        int unclassified = ms.get_unclassified_molecules_count(*p);
        missing[set.path()] = unclassified;
    }

    auto x = std::min_element(missing.begin(), missing.end(), [](const auto &p1, const auto &p2) {
        return p1.second < p2.second;
    });

    return x->first;
}