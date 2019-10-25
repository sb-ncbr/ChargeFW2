//
// Created by krab1k on 24.1.19.
//

#include <vector>
#include <fmt/format.h>
#include <gemmi/pdb.hpp>

#include "chargefw2.h"
#include "../config.h"
#include "pdb.h"
#include "common.h"
#include "bonds.h"
#include "../periodic_table.h"


MoleculeSet PDB::read_file(const std::string &filename) {
    gemmi::Structure structure;
    try {
        structure = gemmi::read_pdb_file(filename);
    }
    catch (std::exception &){
        fmt::print(stderr, "Cannot load structure from file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    auto molecules = std::make_unique<std::vector<Molecule> >();
    auto atoms = std::make_unique<std::vector<Atom>>();

    /* Read first model only */
    auto model = structure.models[0];
    size_t idx = 0;
    for (const auto &chain: model.chains) {
        for (const auto &residue: chain.residues) {
            bool hetatm = residue.het_flag == 'H';
            for (const auto &atom: residue.atoms) {
                double x = atom.pos.x;
                double y = atom.pos.y;
                double z = atom.pos.z;
                int residue_id = residue.seqid.num.value;
                auto element = PeriodicTable::pte().get_element_by_symbol(get_element_symbol(atom.element.name()));

                if (not atom.has_altloc() or not is_already_loaded(*atoms, atom.name, residue_id)) {
                    if ((not hetatm) or
                        (config::read_hetatm and residue.name != "HOH") or
                        (config::read_hetatm and not config::ignore_water)) {
                        atoms->emplace_back(idx, element, x, y, z, atom.name, residue_id, residue.name, chain.name, hetatm);
                        atoms->back()._set_formal_charge(atom.charge);
                        idx++;
                    }
                }
            }
        }
    }

    auto bonds = get_bonds(atoms);
    molecules->emplace_back(sanitize_name(structure.name), std::move(atoms), std::move(bonds));

    return MoleculeSet(std::move(molecules));
}
