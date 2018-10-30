//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <QString>
#include <utility>

class Element {
    int Z_{};
    QString symbol_;
    QString name_;
    float electronegativity_{};
public:
    Element() = default;

    Element(int Z, QString symbol, QString name, float electronegativity) : Z_(Z), symbol_(std::move(symbol)),
                                                                            name_(std::move(name)),
                                                                            electronegativity_(electronegativity) {}

    const QString &symbol() const { return symbol_; }

    float electronegativity() const { return electronegativity_; }

};
