//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <array>
#include <QString>
#include <iostream>

#include "../Element.h"

class Atom {
    int index_{};
    const Element *element_{};
    std::array<double, 3> pos_{};
public:
    Atom() = default;

    Atom(int index, const Element *element, double x, double y, double z);

    int index() const { return index_; }

    const QString &symbol() const { return element_->symbol(); }

    friend std::ostream &operator<<(std::ostream &str, const Atom &atom);

    bool operator==(const Atom &other) const;
};
