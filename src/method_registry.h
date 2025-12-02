#pragma once

#include <memory>
#include <ranges>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <stdexcept>

#include "method.h"


class MethodRegistry {
public:
    using Factory = std::unique_ptr<Method>(*)();

    static void register_factory(std::string name, Factory f) {
        auto& r = map();
        auto [it, ok] = r.emplace(std::move(name), f);
        if (!ok) throw std::runtime_error("Duplicate method name");
    }

    static std::unique_ptr<Method> make(std::string_view name) {
        auto& r = map();
        if (auto it = r.find(std::string{name}); it != r.end()) return it->second();
        return nullptr;
    }

    static std::vector<std::string> names() {
        auto& r = map();
        std::vector<std::string> k; k.reserve(r.size());
        for (const auto& key : r | std::views::keys) k.push_back(key);
        return k;
    }

private:
    using Map = std::unordered_map<std::string, Factory>;
    static Map& map() { static Map m; return m; }
};
