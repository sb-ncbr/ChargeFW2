#include "formats/SDF.h"
#include "structures/MoleculeSet.h"
#include "PeriodicTable.h"


#include <iostream>


int main() {
    SDF reader;
    MoleculeSet m = reader.read_file("/home/krab1k/Research/ChargeFW2/test/set01.sdf");
    m.info();
    return 0;
}