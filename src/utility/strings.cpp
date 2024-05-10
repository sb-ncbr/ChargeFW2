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
