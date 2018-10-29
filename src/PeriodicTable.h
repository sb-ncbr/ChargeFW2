//
// Created by krab1k on 29/10/18.
//

#pragma once

#include <QVector>
#include <QMap>

#include "Element.h"

class PeriodicTable {
    QVector<Element> elements_;
    QMap<QString, Element> symbol_element_;
public:
    static const PeriodicTable &pte();

    const Element getElement(int Z) const { return elements_[Z - 1]; }

    const Element getElement(const QString &symbol) const { return symbol_element_[symbol]; }

    PeriodicTable();
};


