//
// Created by krab1k on 29/10/18.
//

#include <QFile>
#include <QStringList>
#include <QTextStream>
#include <iostream>

#include "Element.h"
#include "PeriodicTable.h"

PeriodicTable::PeriodicTable() {
    QFile file("../data/pte.csv");
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        exit(EXIT_FAILURE);

    QTextStream in(&file);
    // Read header;
    QString line = in.readLine();
    while (!in.atEnd()) {
        line = in.readLine();
        QStringList cols = line.split(',');
        int index = cols[0].toInt();
        QString name = cols[1];
        QString symbol = cols[2];
        float electronegativity = cols[12].toFloat();
        Element element(index, symbol, name, electronegativity);
        elements_.push_back(element);
        symbol_Z_[symbol] = index;
    }
}

const PeriodicTable &PeriodicTable::pte() {
    static PeriodicTable pte;
    return pte;
}

const Element *PeriodicTable::getElement(const QString &symbol) const {
    if(!symbol_Z_.count(symbol)) {
        std::cerr << "No such element " << symbol.toStdString() << std::endl;
        exit(EXIT_FAILURE);
    }
    return getElement(symbol_Z_.at(symbol) - 1);
}
