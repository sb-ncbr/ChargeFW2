#pragma once

#include <string_view>

inline constexpr std::string_view VERSION = "3.0.0-alpha.1";


enum class ExitCode : int {
    Success         = 0,
    InternalError   = 1,
    ParameterError  = 2,
    FileError       = 3
};

[[nodiscard]] constexpr int to_int(ExitCode code) noexcept {
    return static_cast<int>(code);
}


inline constexpr int LARGE_MOLECULE_ATOM_COUNT = 20000;
