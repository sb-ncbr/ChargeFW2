// Created by krab1k on 01/11/18.
//

#include <string>

#include "Classifier.h"
#include "structures/Molecule.h"
#include "structures/MoleculeSet.h"

#include <iostream>

std::string HBOClassifier::get_type(const Atom &atom) const {
    const Molecule *molecule = atom.molecule();
    int max_order = 0;

    for (const Bond &bond: molecule->bonds()) {
        if (bond.hasAtom(atom) && max_order < bond.order())
            max_order = bond.order();
    }
    return std::to_string(max_order);
}
