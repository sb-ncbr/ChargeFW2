//
// Created by krab1k on 24/10/18.
//

#include <map>
#include <iostream>
#include <string>

using std::map;       using std::string;
using std::cout;      using std::endl;

#include "MoleculeSet.h"

void MoleculeSet::info() {
    map<string, int> counts;
    for (auto m: molecules) {
        for (auto a : m.atoms()) {
            if (counts.count(a.symbol()))
                counts[a.symbol()] += 1;
            else
                counts[a.symbol()] = 1;
        }
    }

    for (const auto &[key, val]: counts) {
        cout << key << ": " << val << endl;

    }
}
