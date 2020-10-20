//
// Created by krab1k on 4.1.19.
//

#include <cmath>
#include <fmt/format.h>

#include "chargefw2.h"
#include "charges.h"
#include "statistics.h"


double  RMSD(const Charges &charges1, const Charges &charges2) {
    if (charges1.names() != charges2.names()) {
        throw std::runtime_error("Trying to compare two different sets of charges");
    }

    double total_rmsd = 0;
    for (const auto &name: charges1.names()) {
        std::vector<double> q1 = charges1[name];
        std::vector<double> q2 = charges2[name];
        size_t n = q1.size();
        double sum = 0;
        for (size_t i = 0; i < n; i++) {
            double d = q1[i] - q2[i];
            sum += d * d;
        }

        total_rmsd += std::sqrt(sum / n);
    }
    return total_rmsd / charges1.names().size();
}
