#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <filesystem>
#include <memory>
#include <dlfcn.h>

#include "structures/molecule_set.h"
#include "formats/reader.h"
#include "method.h"
#include "statistics.h"

namespace py = pybind11;
using namespace pybind11::literals;


struct Data {
    MoleculeSet ms;
    Charges reference_charges;
    std::unique_ptr<Parameters> parameters;


    Data(const std::string &input_file, const std::string &ref_chg_file, const std::string &parameter_file);
};


double evaluate(const Data &data, const std::string &method);


Data::Data(const std::string &input_file, const std::string &ref_chg_file, const std::string &parameter_file)
        : reference_charges(Charges(ref_chg_file)) {

    ms = load_molecule_set(input_file);
    if (ms.molecules().empty()) {
        throw std::runtime_error("No molecules were loaded from the input file");
    }

    ms.fulfill_requirements(
            {RequiredFeatures::BOND_DISTANCES, RequiredFeatures::DISTANCE_TREE, RequiredFeatures::BOND_INFO});

    if (not parameter_file.empty()) {
        parameters = std::make_unique<Parameters>(parameter_file);
        auto unclassified = ms.classify_set_from_parameters(*parameters, false, true);
        if (unclassified) {
            throw std::runtime_error("Selected parameters doesn't cover the whole molecule set");
        }
    }
}


double evaluate(const Data &data, const std::string &method_name) {

    auto handle = dlopen(method_name.c_str(), RTLD_LAZY);

    auto get_method_handle = (Method *(*)()) (dlsym(handle, "get_method"));
    if (!get_method_handle) {
        throw std::runtime_error(dlerror());
    }

    auto method = (*get_method_handle)();

    if (method->has_parameters()) {
        method->set_parameters(data.parameters.get());
    }

    Charges charges;
    for (auto &mol: data.ms.molecules()) {
        auto results = method->calculate_charges(mol);
        if (std::any_of(results.begin(), results.end(), [](double chg) { return not isfinite(chg); })) {
            return -1;
        }
        charges.insert(mol.name(), results);
    }

    dlclose(handle);

    return RMSD(data.reference_charges, charges);
}


PYBIND11_MODULE(chargefw2_python, m) {
    m.doc() = "Python binding to ChargeFW2";
    py::class_<Data>(m, "Data")
            .def(py::init<const std::string &, const std::string &, const std::string &>());

    m.def("evaluate", &evaluate, "data"_a, "method"_a, "Evaluate method against reference");
}