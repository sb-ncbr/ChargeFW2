//
// Created by krab1k on 7.1.19.
//

#include <vector>
#include <random>
#include <string>
#include <algorithm>

#include "utility.h"


std::vector<double> generate_random_vector(size_t n, double min, double max) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<double> values(n);
    std::uniform_real_distribution<> dis(min, max);
    std::generate(values.begin(), values.end(), [&]() { return dis(gen); });
    return values;
}


std::string to_lowercase(const std::string &from) {
    std::string result(from);
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}
