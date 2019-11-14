//
// Created by krab1k on 2019-11-14.
//


#include <string>
#include <algorithm>

#include "strings.h"


std::string to_lowercase(const std::string &from) {
    std::string result(from);
    std::transform(result.begin(), result.end(), result.begin(), tolower);
    return result;
}


std::string to_uppercase(const std::string &from) {
    std::string result(from);
    std::transform(result.begin(), result.end(), result.begin(), toupper);
    return result;
}


bool starts_with(const std::string &text, const std::string &prefix) {
    return text.rfind(prefix, 0) == 0;
}


std::string trim(const std::string &text) {
    auto start = text.begin();
    while (start != text.end() and isspace(*start)) {
        start++;
    }
    auto end = text.end();
    do {
        end--;
    } while (std::distance(start, end) > 0 and isspace(*end));

    return std::string(start, end + 1);
}
