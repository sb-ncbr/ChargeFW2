//
// Created by krab1k on 20.11.18.
//

#include <vector>
#include <iostream>
#include <iomanip>
#include <set>
#include <nlopt.hpp>

#include "parameterization.h"
#include "parameters.h"
#include "charges.h"
#include "statistics.h"
#include "utility/utility.h"


double f(const std::vector<double> &x, std::vector<double> &, void *data) {
    auto parameterization = reinterpret_cast<Parameterization*>(data);

    parameterization->method()->parameters()->set_from_vector(x);

    auto charges = Charges();

    for (const auto &molecule: parameterization->set().molecules()) {
        charges.insert(molecule.name(), parameterization->method()->calculate_charges(molecule));
    }

    return RMSD(parameterization->reference_charges(), charges);
}


void Parameterization::parametrize() {

    auto values = method_->parameters()->get_vector();
    size_t n = values.size();

    nlopt::opt opt(nlopt::GN_DIRECT_L, n);
    nlopt::opt local_opt(nlopt::LN_NEWUOA, n);

    std::vector<double> lb(n, 0);
    std::vector<double> ub(n, 3);
    auto x = generate_random_vector(n, 0, 3);

    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    opt.set_min_objective(f, reinterpret_cast<void *>(this));

    opt.set_local_optimizer(local_opt);

    double minf;
    try {
        auto r = opt.optimize(x, minf);
        std::cout << std::setprecision(2) << "found minimum = " << minf << std::endl;
        std::cout << r << std::endl;
    }
    catch (std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }

    std::cout << x << std::endl;
    parameters_->print();
}

Parameterization::Parameterization(const MoleculeSet &ms, boost::shared_ptr<Method> method,
                                   const Charges &reference_charges) : set_{ms}, method_{method},
                                                                       reference_charges_{reference_charges} {

    parameters_ = std::unique_ptr<Parameters>();
    parameters_ = std::make_unique<Parameters>(ms, method);

    method_->set_parameters(parameters_.get());
}
