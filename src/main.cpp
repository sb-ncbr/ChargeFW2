#include <iostream>
#include <fstream>
#include <vector>

#include "formats/SDF.h"
#include "structures/MoleculeSet.h"
#include "PeriodicTable.h"
#include "Parameters.h"
#include "methods/EEM.h"
#include "utility/Utility.h"

int main() {
    SDF reader;
    //MoleculeSet m = reader.read_file("/home/krab1k/Research/Data/NEEMP_Publication/Structures/ideal_valid_all.sdf");
    MoleculeSet m = reader.read_file("/home/krab1k/Research/ChargeFW2/test/set01.sdf");
    const auto p = Parameters("/home/krab1k/Research/ChargeFW2/test/eem.json");
    std::cout << "Parameters:" << std::endl;
    p.print();
    m.classify_atoms_from_parameters(p);
    std::cout << "Set info:" << std::endl;
    m.info();

    auto eem = EEM(&p);

    std::ofstream f("/home/krab1k/Research/ChargeFW2/test/out");

    std::vector<std::vector<double>> grr;

    for(auto &mol: m.molecules()) {
        grr.emplace_back(eem.calculate_charges(mol));
    }

    for(const auto &v: grr) {
        f << v << std::endl;
    }

    f.close();
    return 0;
}