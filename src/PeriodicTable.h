//
// Created by krab1k on 29/10/18.
//

#pragma once

#include <QVector>
#include <QMap>

#include "Element.h"

class PeriodicTable {
    QVector<Element> elements_;
    QMap<QString, int> symbol_Z_;
public:
    static const PeriodicTable &pte();

    const Element &getElement(int Z) const { return elements_[Z - 1]; }

    const Element &getElement(const QString &symbol) const { return elements_[symbol_Z_[symbol] - 1]; }

    PeriodicTable();
};