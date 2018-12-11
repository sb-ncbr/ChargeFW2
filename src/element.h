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
    float covalent_radius_{};
    float vdw_radius_{};
    int group_{};
public:
    Element() = default;

    Element(int Z, std::string symbol, std::string name, float electronegativity, float covalent_radius_,
            float vdw_radius_, int group_)
            : Z_(Z), symbol_(std::move(symbol)), name_(std::move(name)), electronegativity_(electronegativity),
              covalent_radius_(covalent_radius_), vdw_radius_(vdw_radius_), group_(group_) {}

    const std::string &symbol() const { return symbol_; }

    float electronegativity() const { return electronegativity_; }

    float covalent_radius() const { return covalent_radius_; }

    float vdw_radius() const { return vdw_radius_; }

    int group() const { return group_; }

    int valence_electron_count() const;

    int Z() const { return Z_; }
};


