#include <fmt/format.h>
#include <fmt/ostream.h>
#include <nlohmann/json.hpp>
#include <memory>
#include <filesystem>
#include <cstdio>
#include <sys/resource.h>
#include <ctime>
#include <chrono>
#include <unistd.h>
#include <algorithm>

#include "chargefw2.h"
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
#include "candidates.h"
#include "config.h"
#include "parameterization.h"
#include "utility/utility.h"


int main(int argc, char **argv) {
    auto parsed = parse_args(argc, argv);
    check_common_args();

    auto start = std::chrono::system_clock::now();

    auto ext = std::filesystem::path(config::input_file).extension();

    bool is_protein_structure = false;
    std::unique_ptr<Reader> reader;
    if (ext == ".sdf") {
        reader = std::make_unique<SDF>();
    } else if (ext == ".mol2") {
        reader = std::make_unique<Mol2>();
    } else if (ext == ".pdb" or ext == ".ent") {
        reader = std::make_unique<PDB>();
        is_protein_structure = true;
    } else if (ext == ".cif") {
        reader = std::make_unique<mmCIF>();
        is_protein_structure = true;
    } else {
        fmt::print(stderr, "Filetype {} not supported\n", ext);
        exit(EXIT_FILE_ERROR);
    }

    MoleculeSet m = reader->read_file(config::input_file);

    if (config::mode == "info") {
        m.classify_atoms(AtomClassifier::HBO);
        m.info();

    } else if (config::mode == "charges") {
        std::shared_ptr<Method> method = load_method(config::method_name);

        setup_method_options(method, parsed);

        auto p = std::unique_ptr<Parameters>();

        if (method->has_parameters()) {
            std::string par_name;
            if (config::par_file.empty()) {
                par_name = best_parameters(m, method, is_protein_structure, config::permissive_types);
                if (par_name.empty()) {
                    fmt::print(stderr, "No parameters found \n");
                    exit(EXIT_PARAMETER_ERROR);
                }
                fmt::print("Best parameters found: {}\n", par_name);
            } else {
                par_name = config::par_file;
            }

            p = std::make_unique<Parameters>(par_name);
            p->print();

            size_t unclassified = m.classify_set_from_parameters(*p, true, config::permissive_types);
            fmt::print("\nNumber of unclassified molecules: {}\n\n", unclassified);

            method->set_parameters(p.get());
        } else {
            m.classify_atoms(AtomClassifier::PLAIN);
        }

        m.info();
        fmt::print("\n");

        m.fulfill_requirements(method->get_requirements());

        auto charges = Charges();

        charges.set_method_name(config::method_name);

        for (auto &mol: m.molecules()) {
            charges.insert(mol.name(), method->calculate_charges(mol));
        }

        auto txt = TXT();
        std::filesystem::path dir(config::chg_out_dir);
        std::filesystem::path file(config::input_file);
        auto txt_str = file.filename().string() + ".txt";
        txt.save_charges(m, charges, dir / std::filesystem::path(txt_str));

        if (is_protein_structure) {
            auto pqr = PQR();
            auto pqr_str = file.filename().string() + ".pqr";
            pqr.save_charges(m, charges, dir / std::filesystem::path(pqr_str));
        } else {
            auto mol2 = Mol2();
            auto mol2_str = file.filename().string() + ".mol2";
            mol2.save_charges(m, charges, dir / std::filesystem::path(mol2_str));
        }

        if (not config::log_file.empty()) {
            struct rusage usage = {};
            getrusage(RUSAGE_SELF, &usage);

            double utime = usage.ru_utime.tv_sec + static_cast<double>(usage.ru_utime.tv_usec) / 10e6;
            double stime = usage.ru_stime.tv_sec + static_cast<double>(usage.ru_stime.tv_usec) / 10e6;
            double mem = static_cast<double>(usage.ru_maxrss) / 1024;
            auto now = time(nullptr);
            char current_time[100];
            strftime(current_time, 100, "%c", localtime(&now));
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> walltime = end - start;

            auto pid = getpid();
            std::string parameters;
            if (method->parameters()) {
                parameters = method->parameters()->name();
            } else {
                parameters = std::string("None");
            }

            auto log_file = std::fopen(config::log_file.c_str(), "a");

            fmt::print(log_file, "{} [{}]; File: {}; Processed molecules: {}; Method: {}; Parameters: {}\n",
                       current_time, pid, config::input_file, m.molecules().size(), method->name(), parameters);

            fmt::print(log_file,
                       "{} [{}]; Walltime: {:.2f} s; User time: {:.2f} s; System time: {:.2f} s; Peak memory: {:.1f} MB \n",
                       current_time, pid, walltime.count(), utime, stime, mem);
        }
    } else if (config::mode == "best-parameters") {
        std::shared_ptr<Method> method = load_method(config::method_name);

        if (!method->has_parameters()) {
            fmt::print(stderr, "Method uses no parameters\n");
            exit(EXIT_PARAMETER_ERROR);
        }
        auto best = best_parameters(m, method, is_protein_structure, config::permissive_types);
        if (best.empty()) {
            fmt::print("There are no best parameters\n");
        } else {
            fmt::print("Best parameters are: {}\n", best);
        }
    } else if (config::mode == "parameters") {
        m.classify_atoms(AtomClassifier::PLAIN);

        std::shared_ptr<Method> method = load_method(config::method_name);

        setup_method_options(method, parsed);

        Charges reference_charges(config::ref_chg_file);

        auto p = Parameterization(m, method, reference_charges, config::chg_out_dir, config::par_file);
        p.parametrize();

    } else if (config::mode == "suitable-methods") {
        get_suitable_methods(m, is_protein_structure, config::permissive_types);
    } else {
        fmt::print(stderr, "Unknown mode {}\n", config::mode);
        exit(EXIT_PARAMETER_ERROR);
    }

    return 0;
}
