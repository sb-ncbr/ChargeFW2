#include <iostream>
#include <boost/program_options.hpp>
#include <boost/dll/import.hpp>
#include <memory>
#include <ctime>
#include <iomanip>

#include "formats/sdf.h"
#include "structures/molecule_set.h"
#include "parameters.h"
#include "classifier.h"
#include "charges.h"
#include "method.h"
#include "config.h"

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
    try {
        po::store(po::command_line_parser(argc, argv).options(desc).style(
                po::command_line_style::unix_style ^ po::command_line_style::allow_guessing).run(), vm);

        if (vm.count("help")) {
            std::cout << "ChargeFW2 (version " << VERSION << ")" << std::endl;
            std::cout << "by Tomáš Raček (2018)\n" << std::endl;
            std::cout << desc << std::endl;
            exit(EXIT_SUCCESS);
        }
        po::notify(vm);
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    auto mode = vm["mode"].as<std::string>();
    if (mode == "info") {
        if (!vm.count("sdf-file")) {
            std::cerr << "SDF must be provided" << std::endl;
            exit(EXIT_FAILURE);
        }

        auto sdf_name = vm["sdf-file"].as<std::string>();
        SDF reader;
        MoleculeSet m = reader.read_file(sdf_name);

        auto hbo = HBOClassifier();
        m.classify_atoms(hbo);
        m.info();

    } else if (mode == "charges") {
        if (!vm.count("sdf-file")) {
            std::cerr << "SDF must be provided" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (!vm.count("chg-file")) {
            std::cerr << "File where to store charges must be provided" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (!vm.count("method")) {
            std::cerr << "No method selected" << std::endl;
            exit(EXIT_FAILURE);
        }

        auto sdf_name = vm["sdf-file"].as<std::string>();
        auto chg_name = vm["chg-file"].as<std::string>();
        auto method_name = vm["method"].as<std::string>();

        SDF reader;
        MoleculeSet m = reader.read_file(sdf_name);

        boost::shared_ptr<Method> method;

        try {
            method = boost::dll::import<Method>(std::string(INSTALL_DIR) + "/lib/" + method_name, "method",
                                                boost::dll::load_mode::append_decorations);
        } catch (std::exception &e) {
            std::cerr << "Unable to load method " << method_name << std::endl;
            exit(EXIT_FAILURE);
        }

        auto p = std::unique_ptr<Parameters>();

        if (method->has_parameters()) {
            if (!vm.count("par-file")) {
                std::cerr << "File with parameters must be provided" << std::endl;
                exit(EXIT_FAILURE);
            }
            auto par_name = vm["par-file"].as<std::string>();

            p = std::make_unique<Parameters>(par_name);

            std::cout << "Parameters:" << std::endl;
            p->print();

            m.classify_atoms_from_parameters(*p);

            method->set_parameters(p.get());
        } else {
            auto plain = PlainClassifier();
            m.classify_atoms(plain);
        }

        std::cout << "\nSet info:" << std::endl;
        m.info();

        auto charges = Charges();

        clock_t begin = clock();

        for (auto &mol: m.molecules()) {
            charges.insert(mol.name(), method->calculate_charges(mol));
        }

        clock_t end = clock();

        std::cout << "Computation took " << std::setprecision(2) << double(end - begin) / CLOCKS_PER_SEC << " seconds"
                  << std::endl;

        charges.save_to_file(chg_name);

    } else {
        std::cerr << "Unknown mode: " << mode << std::endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}