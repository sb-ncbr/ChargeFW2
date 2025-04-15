#include <cstdio>
#include <dlfcn.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <filesystem>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <tuple>
#include <nlohmann/json.hpp>
#include <string>

#include "charges.h"
#include "exceptions/file_exception.h"
#include "formats/cif.h"
#include "formats/mol2.h"
#include "formats/pqr.h"
#include "formats/txt.h"
#include "method.h"
#include "parameters.h"
#include "structures/molecule_set.h"
#include "formats/reader.h"
#include "config.h"
#include "candidates.h"
#include "utility/strings.h"
#include "utility/install.h"


namespace py = pybind11;
using namespace pybind11::literals;

struct PythonMethodMetadata;


std::map<std::string, std::vector<double>>
calculate_charges(struct Molecules &molecules, const std::string &method_name, std::optional<const std::string> &parameters_name, std::optional<const std::string> &chg_out_dir);

void save_charges(const Molecules &molecules, const Charges &charges, const std::string &filename);

std::vector<PythonMethodMetadata> get_available_methods_python();

std::vector<ParametersMetadata> get_available_parameters(const std::string &method_name);

std::optional<ParametersMetadata> get_best_parameters(struct Molecules &molecules, const std::string &method_name, bool permissive_types);

std::vector<std::tuple<PythonMethodMetadata, std::vector<ParametersMetadata>>> get_suitable_methods_python(struct Molecules &molecules);

py::dict atom_type_count_to_dict(const MoleculeSetStats::AtomTypeCount &atom_type_count);

py::dict molecule_info_to_dict(const MoleculeSetStats &stats);

struct PythonMethodMetadata : public MethodMetadata {
    bool has_parameters;

    PythonMethodMetadata(const MethodMetadata &metadata, bool has_parameters) : MethodMetadata(metadata), has_parameters(has_parameters) {}
};

struct Molecules {
    MoleculeSet ms;

    Molecules(const std::string &filename, bool read_hetatm, bool ignore_water, bool permissive_types);

    std::string input_file;

    [[nodiscard]] size_t length() const;
    [[nodiscard]] MoleculeSetStats info();
};

Molecules::Molecules(const std::string &filename, bool read_hetatm = true, bool ignore_water = true, bool permissive_types = false) {
    config::read_hetatm = read_hetatm;
    config::ignore_water = ignore_water;
    config::permissive_types = permissive_types;
    input_file = filename;

    ms = load_molecule_set(filename);
    if (ms.molecules().empty()) {
        throw std::runtime_error("No molecules were loaded from the input file");
    }
}


size_t Molecules::length() const {
    return ms.molecules().size();
}


MoleculeSetStats Molecules::info() {
    ms.classify_atoms(AtomClassifier::PLAIN);
    return ms.get_stats();
}

py::dict atom_type_count_to_dict(const MoleculeSetStats::AtomTypeCount &atom_type_count) {
        return py::dict(
            py::arg("symbol") = atom_type_count.symbol,
            py::arg("count") = atom_type_count.count
    );
}

py::dict molecule_info_to_dict(const MoleculeSetStats &stats) {
    py::list atom_types_list;
    for (auto &count : stats.atom_type_counts) {
        atom_types_list.append(atom_type_count_to_dict(count));
    }

    return py::dict(
        py::arg("total_molecules") = stats.total_molecules,
        py::arg("total_atoms") = stats.total_atoms,
        py::arg("atom_type_counts") = atom_types_list
    );
}

std::vector<ParametersMetadata> get_available_parameters(const std::string &method_name) {
    std::vector<ParametersMetadata> parameters;
    for (const auto &parameter_file: get_parameter_files()) {
        if (not to_lowercase(parameter_file.filename().string()).starts_with(method_name)) {
            continue;
        }

        auto p = std::make_unique<Parameters>(parameter_file);
        if (method_name == p->method_name()) {
            parameters.emplace_back(p->metadata());
        }
    }
    return parameters;
}

std::vector<std::tuple<PythonMethodMetadata, std::vector<ParametersMetadata>>> get_suitable_methods_python(struct Molecules &molecules) {
    auto suitable = get_suitable_methods(molecules.ms, molecules.ms.has_proteins(), config::permissive_types);
    auto result = std::vector<std::tuple<PythonMethodMetadata, std::vector<ParametersMetadata>>>();
    result.reserve(suitable.size());

    for (auto it = suitable.begin(); it != suitable.end(); ++it) {
        const auto &[method, parameters_list] = *it;

        auto metadata = method->metadata();
        auto has_parameters = method->has_parameters();

        std::vector<ParametersMetadata> parameters_metadata;
        if (has_parameters) {
            for (const auto &parameters : parameters_list) {
                parameters_metadata.emplace_back(parameters->metadata());
            }
        }

        result.emplace_back(PythonMethodMetadata(metadata, has_parameters), parameters_metadata);

    }

    return result;
}

std::vector<PythonMethodMetadata> get_available_methods_python() {
    auto methods = get_available_methods();
    std::vector<PythonMethodMetadata> result;
    for (const auto &method : methods) {
        result.emplace_back(PythonMethodMetadata(method->metadata(), method->has_parameters()));
    }
    return result;
}

std::optional<ParametersMetadata> get_best_parameters(struct Molecules &molecules, const std::string &method_name, bool permissive_types) {
    auto method = load_method(method_name);
    auto parameters = best_parameters(molecules.ms, method, molecules.ms.has_proteins(), permissive_types);

    if (not parameters.has_value()) {
        return std::nullopt;
    }

    return parameters.value()->metadata();
}

std::map<std::string, std::vector<double>>
calculate_charges(struct Molecules &molecules, const std::string &method_name, std::optional<const std::string> &parameters_name, std::optional<const std::string> &chg_out_dir) {
    config::chg_out_dir = chg_out_dir.value_or(".");

    Method* method;
    try {
        method = load_method(method_name);
    } catch (FileException &e) {
        throw std::runtime_error(fmt::format("Failed to load method {}: {}", method_name, e.what()));
    }

    molecules.ms.fulfill_requirements(method->get_requirements());

    std::unique_ptr<Parameters> parameters;
    if (method->has_parameters()) {
        if (not parameters_name.has_value()) {
            throw std::runtime_error(std::string("Method ") + method_name + std::string(" requires parameters"));
        }

        std::string parameter_file = InstallPaths::parametersdir() / (parameters_name.value() + ".json");
        if (not parameter_file.empty()) {
            parameters = std::make_unique<Parameters>(parameter_file);
            auto unclassified = molecules.ms.classify_set_from_parameters(*parameters, false, config::permissive_types);
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

    auto to_save = Charges();
    std::map<std::string, std::vector<double>> charges;
    for (auto &mol: molecules.ms.molecules()) {
        auto results = method->calculate_charges(mol);
        if (std::ranges::any_of(results, [](double chg) noexcept { return not isfinite(chg); })) {
            fmt::print("Incorrect values encountered for: {}. Skipping molecule.\n", mol.name());
        } else {
            to_save.insert(mol.name(), results);
            charges[mol.name()] = results;
        }
    }

    save_charges(molecules, to_save, molecules.input_file);

    return charges;
}

void save_charges(const Molecules &molecules, const Charges &charges, const std::string &filename) {
    std::filesystem::path dir(config::chg_out_dir);
    auto file_path = std::filesystem::path(filename);
    auto ext = file_path.extension().string();

    config::input_file = filename;
    CIF().save_charges(molecules.ms, charges, molecules.input_file);

    auto txt_str = file_path.filename().string() + ".txt";
    TXT().save_charges(molecules.ms, charges, dir / std::filesystem::path(txt_str));

    if (molecules.ms.has_proteins()) {
        auto pqr_str = file_path.filename().string() + ".pqr";
        PQR().save_charges(molecules.ms, charges, dir / std::filesystem::path(pqr_str));
    } else {
        auto mol2_str = file_path.filename().string() + ".mol2";
        Mol2().save_charges(molecules.ms, charges, dir / std::filesystem::path(mol2_str));
    }
}


PYBIND11_MODULE(chargefw2, m) {
    m.doc() = "Python bindings to ChargeFW2";
    py::class_<MoleculeSetStats::AtomTypeCount>(m, "AtomTypeCount")
        .def(py::init<>())
        .def_readwrite("symbol", &MoleculeSetStats::AtomTypeCount::symbol)
        .def_readwrite("count", &MoleculeSetStats::AtomTypeCount::count)
        .def("to_dict", [](const MoleculeSetStats::AtomTypeCount &self) {
            return atom_type_count_to_dict(self);
        });

    py::class_<MoleculeSetStats>(m, "MoleculeSetStats")
        .def(py::init<>())
        .def_readwrite("total_molecules", &MoleculeSetStats::total_molecules)
        .def_readwrite("total_atoms", &MoleculeSetStats::total_atoms)
        .def_readwrite("atom_type_counts", &MoleculeSetStats::atom_type_counts)
        .def("to_dict", [](const MoleculeSetStats &self) {
            return molecule_info_to_dict(self);
        });

    py::class_<Molecules>(m, "Molecules")
        .def(py::init<const std::string &, bool, bool, bool>(), py::arg("input_file"), py::arg("read_hetatm") = true,
                py::arg("ignore_water") = false, py::arg("permissive_types") = false)
        .def("__len__", &Molecules::length)
        .def("info", &Molecules::info);

    py::class_<PythonMethodMetadata>(m, "MethodMetadata")
        .def_readwrite("internal_name", &PythonMethodMetadata::internal_name)
        .def_readwrite("full_name", &PythonMethodMetadata::full_name)
        .def_readwrite("publication", &PythonMethodMetadata::publication)
        .def_readwrite("type", &PythonMethodMetadata::type)
        .def_readwrite("priority", &PythonMethodMetadata::priority)
        .def_readwrite("has_parameters", &PythonMethodMetadata::has_parameters);

    py::class_<ParametersMetadata>(m, "ParametersMetadata")
        .def(py::init<>())
        .def_readwrite("internal_name", &ParametersMetadata::internal_name)
        .def_readwrite("full_name", &ParametersMetadata::full_name)
        .def_readwrite("method", &ParametersMetadata::method)
        .def_readwrite("publication", &ParametersMetadata::publication);

    m.def("get_available_methods", &get_available_methods_python, "Return the list of all available methods");
    m.def("get_available_parameters", &get_available_parameters, "method_name"_a,
          "Return the list of all parameters of a given method");
    m.def("get_best_parameters", &get_best_parameters, "molecules"_a, "method_name"_a, "permissive_types"_a = false,
          "Return the best parameters for a given set of molecules and method name");
    m.def("get_suitable_methods", &get_suitable_methods_python, "molecules"_a, "Get methods and parameters that are suitable for a given set of molecules");
    m.def("calculate_charges", &calculate_charges, "molecules"_a, "method_name"_a, py::arg("parameters_name") = py::none(), py::arg("chg_out_dir") = py::none(),
          "Calculate partial atomic charges for a given molecules and method", py::call_guard<py::gil_scoped_release>());
}
