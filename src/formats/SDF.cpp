//
// Created by krab1k on 24/10/18.
//

#include <QFile>
#include <QString>
#include <QVector>
#include <QTextStream>
#include <iostream>


#include "../structures/Atom.h"
#include "../structures/Bond.h"
#include "../structures/Molecule.h"
#include "../PeriodicTable.h"
#include "SDF.h"

using std::cerr;    using std::endl;

MoleculeSet SDF::read_file(const QString &filename) {
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        cerr << "Cannot open file: " << filename.toStdString() << endl;
        exit(EXIT_FAILURE);
    }

    QTextStream in(&file);

    QString line = in.readLine();

    PeriodicTable pte = PeriodicTable::pte();

    QVector<Molecule> molecules;
    while (!line.isNull()) {
        QString name = line; // Read name_ of the molecule
        line = in.readLine(); // Line with comments
        line = in.readLine(); // Line with comments

        // Read line with counts
        line = in.readLine();

        int n_atoms = line.mid(0, 3).toInt();
        int n_bonds = line.mid(3, 3).toInt();

        QVector<Atom> atoms;
        for (int i = 0; i < n_atoms; i++) {
            line = in.readLine();
            double x = line.mid(0, 10).toDouble();
            double y = line.mid(10, 10).toDouble();
            double z = line.mid(20, 10).toDouble();

            auto atom = Atom(i, &pte.getElement(line.mid(31, 3).trimmed()), x, y, z);

            atoms.push_back(atom);
        }

        QVector<Bond> bonds;
        for (int i = 0; i < n_bonds; i++) {
            line = in.readLine();
            int first = line.mid(0, 3).toInt();
            int second = line.mid(3, 3).toInt();
            int order = line.mid(6, 3).toInt();

            auto bond = Bond(atoms[first - 1], atoms[second - 1], order);
            bonds.push_back(bond);
        }

        auto molecule = Molecule(name, atoms, bonds);
        molecules.push_back(molecule);

        do {
            line = in.readLine();
        } while (line != "$$$$");
        line = in.readLine();

    }
    return MoleculeSet(molecules);
}
