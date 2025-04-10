#pragma once

#include <filesystem>

class InstallPaths {
public:
    static const std::filesystem::path& prefix();
    static std::filesystem::path libdir();
    static std::filesystem::path datadir();
};
