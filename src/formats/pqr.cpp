#include <string>

#include "chargefw2.h"
#include "pqr.h"
#include "../structures/molecule_set.h"
#include "../charges.h"
#include "../exceptions/file_exception.h"


void PQR::save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    auto file = std::fopen(filename.c_str(), "w");
    if (!file) {
        throw FileException(fmt::format("Cannot open file: {}", filename));
    }

    const auto &molecule = ms.molecules()[0];

    try {
        auto chg = charges[molecule.name()];
        for (size_t i = 0; i < molecule.atoms().size(); i++) {
            const auto &atom = molecule.atoms()[i];

            /* Match PDB format weirdness */
            std::string name;
            if (atom.element().symbol().length() == 1 and atom.name().length() < 4) {
                name = " ";
            }
            name.append(atom.name());

            fmt::print(file, "{:<6}{:>5d} {:<4s} {:>3s} {:1s} {:>3d}    {:>8.3f}{:>8.3f}{:>8.3f} {:>6.3f} {:>6.3f}\n",
                       atom.hetatm() ? "HETATM": "ATOM", i + 1, name, atom.residue(), atom.chain_id(), atom.residue_id(),
                       atom.pos()[0], atom.pos()[1], atom.pos()[2], chg[i], atom.element().vdw_radius());
        }
    }
    catch (std::out_of_range &) {
        /* Do nothing */
    }

    fclose(file);
}
