#include <iostream>
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>

#include "formats/SDF.h"
#include "structures/MoleculeSet.h"
#include "PeriodicTable.h"
#include "Parameters.h"
#include "Classifier.h"
#include "Charges.h"
#include "methods/EEM.h"
#include "utility/Utility.h"

int main(int argc, char **argv) {
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "Prints this help")
            ("mode", po::value<std::string>()->required(), "Mode")
            ("sdf-file", po::value<std::string>(), "Input SDF file")
            ("par-file", po::value<std::string>(), "File with parameters (json)")
            ("chg-file", po::value<std::string>(), "File to output charges to");


    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(desc).style(
                po::command_line_style::unix_style ^ po::command_line_style::allow_guessing).run(), vm);

        if (vm.count("help")) {
            std::cout << "ChargeFW2" << std::endl;
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
        if (!vm.count("par-file")) {
            std::cerr << "File with parameters must be provided" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (!vm.count("chg-file")) {
            std::cerr << "File where to store charges must be provided" << std::endl;
            exit(EXIT_FAILURE);
        }

        auto sdf_name = vm["sdf-file"].as<std::string>();
        auto par_name = vm["par-file"].as<std::string>();
        auto chg_name = vm["chg-file"].as<std::string>();

        SDF reader;
        MoleculeSet m = reader.read_file(sdf_name);
        const auto p = Parameters(par_name);

        std::cout << "Parameters:" << std::endl;
        p.print();
        m.classify_atoms_from_parameters(p);
        std::cout << "Set info:" << std::endl;
        m.info();

        auto eem = EEM(&p);
        auto charges = Charges();

        for (auto &mol: m.molecules()) {
            charges.insert(mol.name(), eem.calculate_charges(mol));
        }

        charges.save_to_file(chg_name);

    } else {
        std::cerr << "Unknown mode: " << mode << std::endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}