//
// Created by krab1k on 24/10/18.
//

#include <QMap>
#include <iostream>

using std::cout;      using std::endl;

#include "MoleculeSet.h"

void MoleculeSet::info() const {
    QMap<QString, int> counts;
    for (auto m: molecules_) {
        for (auto a : m.atoms()) {
            counts[a.symbol()] += 1;
        }
    }

    for (auto &key: counts.keys()) {
        cout << key.toStdString() << ": " << counts.value(key) << endl;

    }
}
