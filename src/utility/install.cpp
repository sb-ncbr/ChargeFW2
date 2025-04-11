#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <filesystem>
#include <stdexcept>
#include <climits>

#include "install.h"

namespace {
    std::filesystem::path get_install_prefix() {
        auto install_dir = std::getenv("CHARGEFW2_INSTALL_DIR");
        if (install_dir) {
            return std::filesystem::path(install_dir);
        }

        std::vector<char> buf(PATH_MAX);
        ssize_t len = readlink("/proc/self/exe", buf.data(), buf.size());
        if (len == -1) {
            throw std::runtime_error("Failed to resolve /proc/self/exe");
        }

        return std::filesystem::path(buf.begin(), buf.begin() + len)
               .parent_path().parent_path();
    }
}


const std::filesystem::path& InstallPaths::prefix() {
    static std::filesystem::path cached = get_install_prefix();
    return cached;
}

std::filesystem::path InstallPaths::libdir() {
    return prefix() / "lib";
}


std::filesystem::path InstallPaths::datadir() {
    return prefix() / "share";
}
