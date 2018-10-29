//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <QString>

class Element {
    int Z_;
    QString symbol_;
    QString name_;
public:
    Element() = default;

    Element(int Z, QString symbol, QString name) : Z_(Z), symbol_(symbol), name_(name) {}

    const QString &symbol() const { return symbol_; }
};
