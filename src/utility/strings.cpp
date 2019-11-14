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
