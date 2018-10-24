//
// Created by krab1k on 24/10/18.
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "../structures/Atom.h"
#include "../structures/Bond.h"
#include "../structures/Molecule.h"
#include "SDF.h"

using std::cout;        using std::endl;
using std::getline;     using std::string;
using std::ifstream;    using std::stoi;
using std::stod;        using std::vector;
using std::stringstream;

MoleculeSet SDF::read_file(const std::string &filename) {
    ifstream file(filename);
    if (!file) exit(EXIT_FAILURE);

    string line;

    vector<Molecule> molecules;
    while (getline(file, line)) {
        string name = line; // Read name_ of the molecule
        getline(file, line); // Line with comments
        getline(file, line); // Line with comments

        // Read line with counts
        getline(file, line);

        int n_atoms = stoi(line.substr(0, 3));
        int n_bonds = stoi(line.substr(3, 3));

        vector<Atom> atoms;
        for (int i = 0; i < n_atoms; i++) {
            getline(file, line);
            double x = stod(line.substr(0, 10));
            double y = stod(line.substr(10, 10));
            double z = stod(line.substr(20, 10));

            string symbol;
            // Strips whitespace from symbol
            stringstream(line.substr(31, 3)) >> symbol;

            atoms.emplace_back(Atom(i, symbol, x, y, z));
        }

        vector<Bond> bonds;
        for (int i = 0; i < n_bonds; i++) {
            getline(file, line);
            int first = stoi(line.substr(0, 3));
            int second = stoi(line.substr(0, 3));
            int order = stoi(line.substr(0, 3));

            bonds.emplace_back(Bond(atoms[first - 1], atoms[second - 1], order));
        }

        molecules.emplace_back(Molecule(name, atoms, bonds));

        while (getline(file, line));
    }
    return MoleculeSet(molecules);
}
