//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <array>

#include "../Element.h"

class Atom {
    int index_;
    std::string symbol_;
    std::array<double, 3> pos_;
public:
    Atom(int index, std::string symbol, double x, double y, double z);
    int index() {return index_;}
    std::string symbol() {return symbol_;}
};
