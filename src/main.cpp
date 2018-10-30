#include "formats/SDF.h"
#include "structures/MoleculeSet.h"
#include "PeriodicTable.h"


#include <iostream>

int main() {
    SDF reader;
    MoleculeSet m = reader.read_file("/home/krab1k/Downloads/HOH.sdf");
    m.info();
    auto a = m.molecules()[0].atoms()[1];
    a.print();
    for (auto const &molecule: m.molecules()) {
        for (auto bond: molecule.bonds()) {
            bond.print();
            std::cout << bond.hasAtom(a) << std::endl;
        }
    }
    return 0;
}