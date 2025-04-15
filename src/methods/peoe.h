#pragma once


#include <vector>

#include "../structures/molecule.h"
#include "../method.h"


class PEOE : public Method {
    enum common{dampH};
    enum atom{A, B, C};
public:
    explicit PEOE() : Method({"dampH"}, {"A", "B", "C"}, {},
            {
                {"iters", {"iters", "Number of iterations", "int", "7", {}}}
            }) {}

    [[nodiscard]] const MethodMetadata& metadata() const override;
    
    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;
};
