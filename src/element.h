//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <string>
#include <utility>

class Element {
    int Z_{};
    std::string symbol_;
    std::string name_;
    float electronegativity_{};
public:
    Element() = default;

    Element(int Z, std::string symbol, std::string name, float electronegativity) : Z_(Z), symbol_(std::move(symbol)),
                                                                            name_(std::move(name)),
                                                                            electronegativity_(electronegativity) {}

    const std::string &symbol() const { return symbol_; }

    float electronegativity() const { return electronegativity_; }

};
