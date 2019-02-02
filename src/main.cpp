#include <fmt/format.h>
#include <fmt/ostream.h>
#include <boost/program_options.hpp>
#include <boost/dll/import.hpp>
#include <memory>
#include <ctime>
#include <algorithm>
#include <filesystem>

#include "formats/reader.h"
#include "formats/sdf.h"
#include "formats/mol2.h"
#include "formats/pdb.h"
#include "formats/pqr.h"
#include "formats/mmcif.h"
#include "formats/txt.h"
#include "structures/molecule_set.h"
#include "parameters.h"
#include "charges.h"
#include "method.h"
#include "config.h"
#include "parameterization.h"


std::string best_parameters(MoleculeSet &ms, const boost::shared_ptr<Method> &method);

boost::shared_ptr<Method> load_method(const std::string &method_name);


int main(int argc, char **argv) {
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "Prints this help")
            ("mode", po::value<std::string>()->required(), "Mode")
            ("input-file", po::value<std::string>(), "Input file")
            ("par-file", po::value<std::string>(), "File with parameters (json)")
            ("ref-chg-file", po::value<std::string>(), "File with reference charges")
            ("chg-out-dir", po::value<std::string>(), "Directory to output charges to")
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
            fmt::print("ChargeFW2 (version {})\n", VERSION);
            fmt::print("by Tomáš Raček (2018, 2019)\n");
            fmt::print("{}\n", desc);
            exit(EXIT_SUCCESS);
        }
        po::notify(vm);
    } catch (const std::exception &e) {
        fmt::print(stderr, "{}\n", e.what());
        exit(EXIT_PARAMETER_ERROR);
    }

    if (!vm.count("input-file")) {
        fmt::print(stderr, "Input file must be provided\n");
        exit(EXIT_PARAMETER_ERROR);
    }

    auto input_name = vm["input-file"].as<std::string>();

    auto ext = std::filesystem::path(input_name).extension();

    bool protein_structure = false;
    std::unique_ptr<Reader> reader;
    if (ext == ".sdf") {
        reader = std::make_unique<SDF>();
    } else if (ext == ".mol2") {
        reader = std::make_unique<Mol2>();
    } else if (ext == ".pdb") {
        reader = std::make_unique<PDB>();
        protein_structure = true;
    } else if (ext == ".cif") {
        reader = std::make_unique<mmCIF>();
        protein_structure = true;
    } else {
        fmt::print(stderr, "Filetype {} not supported\n", ext);
        exit(EXIT_FILE_ERROR);
    }

    MoleculeSet m = reader->read_file(input_name);

    auto mode = vm["mode"].as<std::string>();
    if (mode == "info") {
        m.classify_atoms(AtomClassifier::HBO);
        m.info();

    } else if (mode == "charges") {
        if (!vm.count("chg-out-dir")) {
            fmt::print(stderr, "Directory where to store charges must be provided\n");
            exit(EXIT_PARAMETER_ERROR);
        }

        if (!vm.count("method")) {
            fmt::print(stderr, "No method selected.\n");
            exit(EXIT_PARAMETER_ERROR);
        }

        auto chg_out_dir = vm["chg-out-dir"].as<std::string>();
        auto method_name = vm["method"].as<std::string>();

        boost::shared_ptr<Method> method = load_method(method_name);

        po::options_description method_options("Method options");
        for (const auto &[opt, info]: method->get_options()) {
            if (!info.choices.empty()) {
                if (std::find(info.choices.begin(), info.choices.end(), info.default_value) == info.choices.end()) {
                    fmt::print(stderr, "Default value: {} not in possible choices\n", info.default_value);
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
            if (vm.count(opt_name)) {
                std::string val = vm[opt_name].as<std::string>();
                if (!info.choices.empty()) {
                    if (std::find(info.choices.begin(), info.choices.end(), val) == info.choices.end()) {
                        fmt::print(stderr, "Provided value: {} not in possible choices\n", val);
                        exit(EXIT_INTERNAL_ERROR);
                    }
                }
                method->set_option_value(opt, val);
            } else {
                method->set_option_value(opt, info.default_value);
            }
        }

        auto p = std::unique_ptr<Parameters>();

        if (method->has_parameters()) {
            std::string par_name;
            if (!vm.count("par-file")) {
                par_name = best_parameters(m, method);
                if (par_name.empty()) {
                    fmt::print(stderr, "No parameters found \n");
                    exit(EXIT_PARAMETER_ERROR);
                }
                fmt::print("Best parameters found: {}\n", par_name);
            } else {
                par_name = vm["par-file"].as<std::string>();
            }

            p = std::make_unique<Parameters>(par_name);

            fmt::print("Parameters:\n");
            p->print();

            size_t unclassified = m.classify_set_from_parameters(*p);
            fmt::print("\nNumber of unclassified molecules: {}\n\n", unclassified);

            method->set_parameters(p.get());
        } else {
            m.classify_atoms(AtomClassifier::PLAIN);
        }

        m.info();
        fmt::print("\n");

        m.fulfill_requirements(method->get_requirements());

        auto charges = Charges();

        charges.set_method_name(method_name);

        clock_t begin = clock();

        for (auto &mol: m.molecules()) {
            charges.insert(mol.name(), method->calculate_charges(mol));
        }

        clock_t end = clock();

        fmt::print("Computation took {:.2f} seconds\n", double(end - begin) / CLOCKS_PER_SEC);

        auto txt = TXT();
        std::filesystem::path dir(chg_out_dir);
        std::filesystem::path file(input_name);
        auto txt_str = file.filename().string() + ".txt";
        txt.save_charges(m, charges, dir / std::filesystem::path(txt_str));

        if (protein_structure) {
            auto pqr = PQR();
            auto pqr_str = file.filename().string() + ".pqr";
            pqr.save_charges(m, charges, dir / std::filesystem::path(pqr_str));
        } else {
            auto mol2 = Mol2();
            auto mol2_str = file.filename().string() + ".mol2";
            mol2.save_charges(m, charges, dir / std::filesystem::path(mol2_str));
        }


    } else if (mode == "best-parameters") {
        if (!vm.count("method")) {
            fmt::print(stderr, "No method selected\n");
            exit(EXIT_PARAMETER_ERROR);
        }

        auto method_name = vm["method"].as<std::string>();

        boost::shared_ptr<Method> method = load_method(method_name);

        if (!method->has_parameters()) {
            fmt::print(stderr, "Method uses no parameters\n");
            exit(EXIT_PARAMETER_ERROR);
        }

        fmt::print("Best parameters are: {}\n", best_parameters(m, method));

    } else if (mode == "parameters") {

        if (!vm.count("ref-chg-file")) {
            fmt::print(stderr, "File with reference charges must be provided");
            exit(EXIT_FAILURE);
        }

        if (!vm.count("chg-out-dir")) {
            fmt::print(stderr, "Directory where to store charges must be provided");
            exit(EXIT_FAILURE);
        }

        if (!vm.count("par-file")) {
            fmt::print(stderr, "File where to store parameters must be provided");
            exit(EXIT_FAILURE);
        }

        if (!vm.count("method")) {
            fmt::print(stderr, "No method selected");
            exit(EXIT_FAILURE);
        }

        auto chg_out_dir = vm["chg-out-dir"].as<std::string>();
        auto par_name = vm["par-file"].as<std::string>();
        auto ref_charge_name = vm["ref-chg-file"].as<std::string>();
        auto method_name = vm["method"].as<std::string>();

        m.classify_atoms(AtomClassifier::PLAIN);

        boost::shared_ptr<Method> method = load_method(method_name);

        Charges reference_charges(ref_charge_name);

        auto p = Parameterization(m, method, reference_charges, chg_out_dir, par_name);
        p.parametrize();

    } else {
        fmt::print(stderr, "Unknown mode {}\n", mode);
        exit(EXIT_PARAMETER_ERROR);
    }

    return 0;
}


std::string best_parameters(MoleculeSet &ms, const boost::shared_ptr<Method> &method) {
    std::string best_name;
    size_t best_unclassified = ms.molecules().size();

    for (const auto &set: std::filesystem::directory_iterator(std::string(INSTALL_DIR) + "/share/parameters")) {
        auto p = std::make_unique<Parameters>(set.path());

        if (method->name() != p->method_name())
            continue;

        size_t unclassified = ms.classify_set_from_parameters(*p, false);

        // If all molecules are covered by the parameters, we found our best
        if (!unclassified) {
            return set.path();
        }

        if (unclassified < best_unclassified) {
            best_unclassified = unclassified;
            best_name = set.path();
        }
    }

    return best_name;
}


boost::shared_ptr<Method> load_method(const std::string &method_name) {
    try {
        return boost::dll::import<Method>(std::string(INSTALL_DIR) + "/lib/" + method_name, "method",
                                          boost::dll::load_mode::append_decorations);
    } catch (std::exception &) {
        fmt::print(stderr, "Unable to load method {}\n", method_name);
        exit(EXIT_PARAMETER_ERROR);
    }
}
