//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <array>
#include <QString>

#include "../Element.h"

class Atom {
    int index_{};
    Element element_;
    std::array<double, 3> pos_{};
public:
    Atom() = default;

    Atom(int index, Element element, double x, double y, double z);

    int index() { return index_; }

    QString symbol() { return element_.symbol(); }
};
