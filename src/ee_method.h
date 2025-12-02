#pragma once

#include <string>
#include <vector>
#include <map>
#include <functional>
#include <Eigen/Core>

#include "method.h"

class EEMethod : public Method {
    [[nodiscard]] std::map<std::string, MethodOption>
    augment_options(std::map<std::string, MethodOption> options) const {
        options["type"] = {"type", "Type of a solver", "str", "full", {"full", "cutoff", "cover"}};
        options["radius"] = {"radius", "Radius for cutoff", "double", "12", {}};
        return options;
    }

public:
    EEMethod(std::vector<std::string> common, std::vector<std::string> atom,
             std::vector<std::string> bond, std::map<std::string, MethodOption> options) :
        Method(std::move(common), std::move(atom), std::move(bond), augment_options(
                   std::move(options))) {
    }

    [[nodiscard]] bool is_suitable_for_large_molecule() const override;

    [[nodiscard]] Eigen::VectorXd solve_EE(const Molecule& molecule,
                                           const std::function<Eigen::VectorXd(const std::vector<const Atom*>&, double)>
                                           &) const;

    [[nodiscard]] std::vector<RequiredFeatures> get_requirements() const override {
        return {RequiredFeatures::DISTANCE_TREE};
    }
};
