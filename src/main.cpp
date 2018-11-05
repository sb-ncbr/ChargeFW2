#include <iostream>

#include "formats/SDF.h"
#include "structures/MoleculeSet.h"
#include "PeriodicTable.h"
#include "Parameters.h"


int main() {
    SDF reader;
    MoleculeSet m = reader.read_file("/home/krab1k/Research/ChargeFW2/test/set01.sdf");
    auto p = Parameters("/home/krab1k/Research/ChargeFW2/test/mpeoe.json");
    p.print();
    m.classify_atoms_from_parameters(p);
    m.info();

    std::cout << p.atom()["b"](m.molecules()[0].atoms()[0]) << std::endl;
    return 0;
}