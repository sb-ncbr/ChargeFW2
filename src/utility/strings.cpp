#include <string>
#include <algorithm>
#include <ranges>
#include "strings.h"


std::string to_lowercase(const std::string &from) {
    auto to_lower = [](unsigned char ch) noexcept { return std::tolower(ch); };

    auto transformed = from | std::views::transform(to_lower);
    return std::string(transformed.begin(), transformed.end());
}


std::string to_uppercase(const std::string &from) {
    auto to_upper = [](unsigned char ch) noexcept { return std::toupper(ch); } ;

    auto transformed = from | std::views::transform(to_upper);
    return std::string(transformed.begin(), transformed.end());
}


std::string trim(const std::string &text) {
    auto not_space = [](unsigned char c) noexcept { return !std::isspace(c); };

    auto first = std::ranges::find_if(text, not_space);
    auto last = std::ranges::find_if(std::ranges::reverse_view(text), not_space).base();

    if (first >= last) {
        return {};
    }

    return std::string(first, last);
}
