//
// Created by krab1k on 24/10/18.
//

#include <QFile>
#include <QString>
#include <QTextStream>
#include <iostream>
#include <tuple>
#include <vector>
#include <memory>


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

    auto molecules = std::make_unique<std::vector<Molecule> >();
    while (!line.isNull()) {
        QString name = line; // Read name_ of the molecule
        line = in.readLine(); // Line with comments
        line = in.readLine(); // Line with comments

        // Read line with counts
        line = in.readLine();

        size_t n_atoms = line.mid(0, 3).toUInt();
        size_t n_bonds = line.mid(3, 3).toUInt();

        auto atoms = std::make_unique<std::vector<Atom> >();
        atoms->reserve(n_atoms);

        for (unsigned i = 0; i < n_atoms; i++) {
            line = in.readLine();
            double x = line.mid(0, 10).toDouble();
            double y = line.mid(10, 10).toDouble();
            double z = line.mid(20, 10).toDouble();

            auto element = PeriodicTable::pte().getElement(line.mid(31, 3).trimmed());

            atoms->emplace_back(i, element, x, y, z);
        }

        auto bonds = std::make_unique<std::vector<Bond> >();
        bonds->reserve(n_bonds);

        for (unsigned i = 0; i < n_bonds; i++) {
            line = in.readLine();
            int first = line.mid(0, 3).toInt();
            int second = line.mid(3, 3).toInt();
            int order = line.mid(6, 3).toInt();

            bonds->emplace_back(&((*atoms)[first - 1]), &((*atoms)[second - 1]), order);
        }

        molecules->emplace_back(name, std::move(atoms), std::move(bonds));

        do {
            line = in.readLine();
        } while (line != "$$$$");
        line = in.readLine();

    }
    return MoleculeSet(std::move(molecules));
}
