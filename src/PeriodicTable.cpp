//
// Created by krab1k on 29/10/18.
//

#include <QFile>
#include <QStringList>
#include <QTextStream>
#include <QDebug>

#include "Element.h"
#include "PeriodicTable.h"

PeriodicTable::PeriodicTable() {
    QFile file("../data/pte.csv");
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
        exit(EXIT_FAILURE);

    QTextStream in(&file);
    // Read header;
    QString line = in.readLine();
    while(!in.atEnd()) {
        line = in.readLine();
        QStringList cols = line.split(',');
        int index = cols[0].toInt();
        QString name = cols[1];
        QString symbol = cols[2];
        Element element(index, symbol, name);
        elements_.push_back(element);
        symbol_element_[symbol] = element;
    }
}

const PeriodicTable &PeriodicTable::pte() {
    static PeriodicTable pte;
    return pte;
}
