//
// Created by krab1k on 4.1.19.
//

#include <cmath>
#include <fmt/format.h>
#include <numeric>

#include "chargefw2.h"
#include "charges.h"
#include "statistics.h"


double RMSD(const Charges &charges1, const Charges &charges2) {
    if (charges1.names() != charges2.names()) {
        throw std::runtime_error("Trying to compare two different sets of charges");
    }

    double sum = 0;
    size_t n = 0;
    for (const auto &name: charges1.names()) {
        std::vector<double> q1 = charges1[name];
        std::vector<double> q2 = charges2[name];
        n += q1.size();
        for (size_t i = 0; i < q1.size(); i++) {
            double d = q1[i] - q2[i];
            sum += d * d;
        }
    }
    return std::sqrt(sum / n);
}


double Pearson2(const Charges &charges1, const Charges &charges2) {
    if (charges1.names() != charges2.names()) {
        throw std::runtime_error("Trying to compare two different sets of charges");
    }

    double mx = 0;
    double my = 0;
    size_t sizes = 0;
    for (const auto &name: charges1.names()) {
        std::vector<double> q1 = charges1[name];
        std::vector<double> q2 = charges2[name];
        size_t n = q1.size();
        mx += std::accumulate(q1.begin(), q1.end(), 0.0);
        my += std::accumulate(q2.begin(), q2.end(), 0.0);
        sizes += n;
    }

    mx /= sizes;
    my /= sizes;

    double num = 0;
    double den1 = 0;
    double den2 = 0;

    for (const auto &name: charges1.names()) {
        std::vector<double> q1 = charges1[name];
        std::vector<double> q2 = charges2[name];
        size_t n = q1.size();

        for (size_t i = 0; i < n; i++) {
            num += (q1[i] - mx) * (q2[i] - my);
            den1 += (q1[i] - mx) * (q1[i] - mx);
            den2 += (q2[i] - my) * (q2[i] - my);
        }
    }

    return (num * num) / (den1 * den2);
}


double D_max(const Charges &charges1, const Charges &charges2) {
    if (charges1.names() != charges2.names()) {
        throw std::runtime_error("Trying to compare two different sets of charges");
    }
    double dmax = -1;
    for (const auto &name: charges1.names()) {
        std::vector<double> q1 = charges1[name];
        std::vector<double> q2 = charges2[name];
        for (size_t i = 0; i < q1.size(); i++) {
            if (fabs(q1[i] - q2[i]) > dmax) {
                dmax = fabs(q1[i] - q2[i]);
            }
        }
    }
    return dmax;
}


double D_avg(const Charges &charges1, const Charges &charges2) {
    if (charges1.names() != charges2.names()) {
        throw std::runtime_error("Trying to compare two different sets of charges");
    }

    double sum = 0;
    size_t n  = 0;
    for (const auto &name: charges1.names()) {
        std::vector<double> q1 = charges1[name];
        std::vector<double> q2 = charges2[name];
        n += q1.size();
        for (size_t i = 0; i < q1.size(); i++) {
            sum += fabs(q1[i] - q2[i]);

        }
    }
    return sum / n;
}
