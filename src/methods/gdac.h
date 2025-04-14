#pragma once

#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class GDAC : public Method {
    enum atom{A, B};
public:
    explicit GDAC() : Method({}, {"A", "B"}, {},
            {
                {"iters", {"iters", "Number of iterations", "int", "7", {}}}
            }) {}

    [[nodiscard]] const MethodMetadata& metadata() const override;
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
