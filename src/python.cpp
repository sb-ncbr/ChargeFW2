#include <cstdio>
#include <fmt/core.h>
#include <fmt/format.h>
#include <dlfcn.h>
#include <filesystem>
#include <fstream>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <nlohmann/json.hpp>

#include "structures/molecule_set.h"
#include "formats/reader.h"
#include "config.h"
#include "candidates.h"
#include "utility/strings.h"


namespace fs = std::filesystem;
namespace py = pybind11;
using namespace pybind11::literals;


std::map<std::string, std::vector<double>>
calculate_charges(struct Molecules &molecules, const std::string &method_name, std::optional<const std::string> &parameters_name);

std::vector<std::string> get_available_methods();

std::vector<std::string> get_available_parameters(const std::string &method_name);

std::vector<std::tuple<std::string, std::vector<std::string>>> get_sutaible_methods_python(struct Molecules &molecules);

struct MoleculeInfo {
    struct AtomTypeCount {
        std::string symbol;
        int count;

        [[nodiscard]] py::dict to_dict() const;
    };

    size_t total_molecules;
    size_t total_atoms;
    std::vector<AtomTypeCount> atom_type_counts;

    [[nodiscard]] py::dict to_dict() const;
};

py::dict MoleculeInfo::to_dict() const {
    py::list atom_types_list;
    for (const auto &count : atom_type_counts) {
        atom_types_list.append(count.to_dict());
    }

    return py::dict(
        py::arg("total_molecules") = total_molecules,
        py::arg("total_atoms") = total_atoms,
        py::arg("atom_type_counts") = atom_types_list
    );
}

py::dict MoleculeInfo::AtomTypeCount::to_dict() const {
        return py::dict(
            py::arg("symbol") = symbol,
            py::arg("count") = count
    );
}

struct Molecules {
    MoleculeSet ms;

    Molecules(const std::string &filename, bool read_hetatm, bool ignore_water);

    [[nodiscard]] size_t length() const;
    [[nodiscard]] MoleculeInfo info();
};

Molecules::Molecules(const std::string &filename, bool read_hetatm = true, bool ignore_water = true) {
    config::read_hetatm = read_hetatm;
    config::ignore_water = ignore_water;
    ms = load_molecule_set(filename);
    if (ms.molecules().empty()) {
        throw std::runtime_error("No molecules were loaded from the input file");
    }

    ms.fulfill_requirements({RequiredFeatures::DISTANCE_TREE, RequiredFeatures::BOND_DISTANCES});
}


size_t Molecules::length() const {
    return ms.molecules().size();
}


std::vector<std::string> get_available_methods() {
    std::vector<std::string> results;
    std::string filename = fs::path(INSTALL_DIR) / "share" / "methods.json";
    using json = nlohmann::json;
    json j;
    std::ifstream f(filename);
    if (!f) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    f >> j;
    f.close();

    for (const auto &method_info: j["methods"]) {
        auto method_name = method_info["internal_name"].get<std::string>();
        results.emplace_back(method_name);
    }

    return results;
}


MoleculeInfo Molecules::info() {
    ms.classify_atoms(AtomClassifier::PLAIN);

    MoleculeInfo result;
    std::map<size_t, int> counts;
    size_t n_atoms = 0;

    for (const auto &m: ms.molecules()) {
        for (auto &a : m.atoms()) {
            counts[a.type()] += 1;
            n_atoms++;
        }
    }

    result.total_molecules = length();
    result.total_atoms = n_atoms;

    auto atom_types = ms.atom_types();

    if (atom_types.size() > 0) {
        for (auto &[key, val] : counts) {
            auto symbol =  std::get<0>(atom_types[key]);
            
            result.atom_type_counts.push_back({
                .symbol = symbol,
                .count = val
            });
        }
    }
    
    return result;
}

std::vector<std::string> get_available_parameters(const std::string &method_name) {
    std::vector<std::string> parameters;
    for (const auto &parameter_file: get_parameter_files()) {
        if (not to_lowercase(parameter_file.filename().string()).starts_with(method_name)) {
            continue;
        }

        auto p = std::make_unique<Parameters>(parameter_file);
        if (method_name == p->method_name()) {
            parameters.emplace_back(parameter_file.stem().string());
        }
    }
    return parameters;
}

std::vector<std::tuple<std::string, std::vector<std::string>>> get_sutaible_methods_python(struct Molecules &molecules) {
    return get_suitable_methods(molecules.ms, molecules.ms.has_proteins(), false);
}


std::map<std::string, std::vector<double>>
calculate_charges(struct Molecules &molecules, const std::string &method_name, std::optional<const std::string> &parameters_name) {
    std::string method_file = fs::path(INSTALL_DIR) / "lib" / ("lib" + method_name + ".so");
    auto handle = dlopen(method_file.c_str(), RTLD_LAZY);

    auto get_method_handle = reinterpret_cast<Method *(*)()>(dlsym(handle, "get_method"));
    if (!get_method_handle) {
        throw std::runtime_error(dlerror());
    }

    auto method = (*get_method_handle)();
    std::unique_ptr<Parameters> parameters;
    if (method->has_parameters()) {
        if (!parameters_name.has_value()) {
            throw std::runtime_error(std::string("Method ") + method_name + std::string(" requires parameters"));
        }

        std::string parameter_file = fs::path(INSTALL_DIR) / "share" / "parameters" / (parameters_name.value() + ".json");
        if (not parameter_file.empty()) {
            parameters = std::make_unique<Parameters>(parameter_file);
            auto unclassified = molecules.ms.classify_set_from_parameters(*parameters, false, true);
            if (unclassified) {
                throw std::runtime_error("Selected parameters doesn't cover the whole molecule set");
            }
        }
    }

    method->set_parameters(parameters.get());

    /* Use only default values */
    for (const auto &[opt, info]: method->get_options()) {
        method->set_option_value(opt, info.default_value);
    }

    std::map<std::string, std::vector<double>> charges;
    for (auto &mol: molecules.ms.molecules()) {
        auto results = method->calculate_charges(mol);
        if (std::any_of(results.begin(), results.end(), [](double chg) { return not isfinite(chg); })) {
            fmt::print("Incorrect values encoutened for: {}. Skipping molecule.\n", mol.name());
        } else {
            charges[mol.name()] = results;
        }
    }

    dlclose(handle);

    return charges;
}


PYBIND11_MODULE(chargefw2, m) {
    m.doc() = "Python bindings to ChargeFW2";
    py::class_<MoleculeInfo::AtomTypeCount>(m, "AtomTypeCount")
        .def(py::init<>())
        .def_readwrite("symbol", &MoleculeInfo::AtomTypeCount::symbol)
        .def_readwrite("count", &MoleculeInfo::AtomTypeCount::count)
        .def("to_dict", &MoleculeInfo::AtomTypeCount::to_dict);

    py::class_<MoleculeInfo>(m, "MoleculeInfo")
        .def(py::init<>())
        .def_readwrite("total_molecules", &MoleculeInfo::total_molecules)
        .def_readwrite("total_atoms", &MoleculeInfo::total_atoms)
        .def_readwrite("atom_type_counts", &MoleculeInfo::atom_type_counts)
        .def("to_dict", &MoleculeInfo::to_dict);

    py::class_<Molecules>(m, "Molecules")
        .def(py::init<const std::string &, bool, bool>(), py::arg("input_file"), py::arg("read_hetatm") = true,
                py::arg("ignore_water") = false)
        .def("__len__", &Molecules::length)
        .def("info", &Molecules::info);

    m.def("get_available_methods", &get_available_methods, "Return the list of all available methods");
    m.def("get_available_parameters", &get_available_parameters, "method_name"_a,
          "Return the list of all parameters of a given method");
    m.def("get_suitable_methods", &get_sutaible_methods_python, "molecules"_a, "Get methods and parameters that are suitable for a given set of molecules");
    m.def("calculate_charges", &calculate_charges, "molecules"_a, "method_name"_a, py::arg("parameters_name") = py::none(),
          "Calculate partial atomic charges for a given molecules and method", py::call_guard<py::gil_scoped_release>());
}
