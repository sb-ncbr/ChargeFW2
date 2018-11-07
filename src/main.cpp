#include <iostream>
#include <vector>

#include "formats/SDF.h"
#include "structures/MoleculeSet.h"
#include "PeriodicTable.h"
#include "Parameters.h"
#include "methods/EEM.h"
#include "utility/Utility.h"

int main() {
    SDF reader;
    MoleculeSet m = reader.read_file("/home/krab1k/Research/ChargeFW2/test/set01.sdf");
    const auto p = Parameters("/home/krab1k/Research/ChargeFW2/test/eem.json");
    std::cout << "Parameters:" << std::endl;
    p.print();
    m.classify_atoms_from_parameters(p);
    std::cout << "Set info:" << std::endl;
    m.info();
    auto atom = m.molecules()[0].atoms()[0];

    auto eem = EEM(&p);
    for(auto &mol: m.molecules()) {
        std::cout << mol.name() << std::endl;
        auto res = eem.calculate_charges(mol);
        std::cout << res << std::endl;
    }
    return 0;
}