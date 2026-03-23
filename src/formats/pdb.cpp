#include <vector>
#include <print>
#include <format>
#include <gemmi/pdb.hpp>

#include "../config.h"
#include "pdb.h"
#include "common.h"
#include "bonds.h"
#include "../periodic_table.h"
#include "../utility/exceptions.h"


MoleculeSet PDB::read_file(const std::string &filename) {
    gemmi::Structure structure;
    try {
        structure = gemmi::read_pdb_file(filename);
    }
    catch (std::exception &e) {
        throw FileException(std::format("Cannot load structure: {}", e.what()));
    }

    auto molecules = std::make_unique<std::vector<Molecule>>();
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

                const Element *element;
                try {
                    element = PeriodicTable::pte().get_element_by_symbol(get_element_symbol(atom.element.name()));
                } catch (std::exception &e) {
                    std::println(stderr, "Error when reading {}: {}", structure.name, e.what());
                    /* Return empty set */
                    return MoleculeSet(std::move(molecules));
                }
                if (keep_atom(atom, residue)) {
                    atoms->emplace_back(idx, element, x, y, z, atom.name, residue_id, residue.name, chain.name, hetatm);
                    atoms->back()._set_formal_charge(atom.charge);
                    idx++;
                }
            }
        }
    }

    if (atoms->empty()) {
        throw FileException(std::format("Error when reading {}: No atoms were loaded", structure.name));
    }

    auto bonds = get_bonds(atoms);
    std::string name = structure.name;
    auto it = structure.info.find("_entry.id");
    if (it != structure.info.end()) {
        name = it->second;
    }
    molecules->emplace_back(sanitize_name(name), std::move(atoms), std::move(bonds));

    return MoleculeSet(std::move(molecules));
}

PDB::PDB() = default;
