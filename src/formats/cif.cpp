#include <string>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <filesystem>

#include <gemmi/align.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/to_mmcif.hpp>

#include "chargefw2.h"
#include "cif.h"
#include "../config.h"
#include "../utility/strings.h"
#include "../exceptions/file_exception.h"


static std::string convert_bond_order_to_mmcif_value_order_string(int order) {
    // See link for bond order values for PDBx/mmCIF category _chem_comp_bond.value_order
    // https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_bond.value_order.html#papwtenum
    switch (order) {
        case 1:
            return "SING";
        case 2:
            return "DOUB";
        case 3:
            return "TRIP";
        case 4:
            return "AROM";
        default:
            return ".";  // unknown
    }
}

static void append_charges_to_block(const Molecule &molecule, const Charges &charges, gemmi::cif::Block &block) {
    const std::string sb_ncbr_partial_atomic_charges_meta_prefix = "_sb_ncbr_partial_atomic_charges_meta.";
    const std::string sb_ncbr_partial_atomic_charges_prefix = "_sb_ncbr_partial_atomic_charges.";
    
    const std::vector<std::string> sb_ncbr_partial_atomic_charges_meta_attributes = {
        "id",
        "type",
        "method",
    };

    const std::vector<std::string> sb_ncbr_partial_atomic_charges_attributes = {
        "type_id",
        "atom_id",
        "charge",
    };

    const auto& atom_charges = charges[molecule.name()];

    // _sb_ncbr_partial_atomic_charges_meta
    auto& metadata_loop = block.init_loop(sb_ncbr_partial_atomic_charges_meta_prefix, sb_ncbr_partial_atomic_charges_meta_attributes);
    const auto id = "1";
    const auto type = "empirical";
    const auto method = fmt::format("'{}/{}'", charges.method_name(), charges.parameters_name());
    metadata_loop.add_row({
        id,
        type,
        method,
    });

    // _sb_ncbr_partial_atomic_charges  
    auto& charges_loop = block.init_loop(sb_ncbr_partial_atomic_charges_prefix, sb_ncbr_partial_atomic_charges_attributes);
    for (size_t i = 0; i < molecule.atoms().size(); ++i) {
        const auto &atom = molecule.atoms()[i];
        const auto atomId = fmt::format("{}", atom.index() + 1);
        const auto charge = fmt::format("{: 1.4f}", atom_charges[i]);
        charges_loop.add_row({
            id,
            atomId,
            charge,
        });
    }
}

static void filter_out_altloc_atoms(gemmi::cif::Block &block) {
    auto structure = gemmi::make_structure_from_block(block);
    if (structure.models.empty()) {
        throw std::runtime_error("Not enough information to create a structure");
    }

    // filter out the altloc atoms
    // retains the _atom_site.id order
    auto &model = structure.models[0];
    for (gemmi::Chain &chain: model.chains) {
        for (gemmi::Residue &residue: chain.residues) {
            bool hetatm = residue.het_flag == 'H';
            for (auto it = residue.atoms.begin(); it != residue.atoms.end(); ) {
                const auto& atom = *it;
                if (not atom.has_altloc() or atom.altloc == 'A') {
                    if ((not hetatm) or
                        (config::read_hetatm and residue.name != "HOH") or
                        (config::read_hetatm and not config::ignore_water)) {
                        ++it;
                    }
                } else {
                    residue.atoms.erase(it);
                }
            }
        }
    }

    gemmi::update_mmcif_block(structure, block);
}

static void generate_mmcif_from_block(gemmi::cif::Block &block, const MoleculeSet &ms, const Charges &charges) {
    const Molecule &molecule = ms.molecules()[0];
    
    filter_out_altloc_atoms(block);
    append_charges_to_block(molecule, charges, block);
    
    const std::filesystem::path out_dir{config::chg_out_dir};
    const std::string molecule_name = to_lowercase(molecule.name());
    const std::string out_filename = molecule_name + ".fw2.cif";
    std::ofstream out_stream{{out_dir / out_filename}};

    // remove pesky _chem_comp category >:(
    block.find_mmcif_category("_chem_comp.").erase();

    gemmi::cif::write_cif_block_to_stream(out_stream, block);
}

static void generate_mmcif_from_cif_file(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    auto document = gemmi::cif::read_file(filename);
    auto& block = document.sole_block();
    
    generate_mmcif_from_block(block, ms, charges);
}

static void generate_mmcif_from_pdb_file(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    auto structure = gemmi::read_pdb_file(filename);

    if (structure.models.empty() || structure.models[0].chains.empty()) {
        throw FileException("No models or no chains in PDB file.");
    }

    gemmi::setup_entities(structure);
    gemmi::assign_label_seq_id(structure, false);

    auto block = gemmi::make_mmcif_block(structure);

    generate_mmcif_from_block(block, ms, charges);
}

static void generate_mmcif_from_atom_and_bond_data(const MoleculeSet &ms, const Charges &charges) {
    const std::string atom_site_prefix = "_atom_site.";
    const std::string chem_comp_prefix = "_chem_comp.";    
    const std::string chem_comp_bond_prefix = "_chem_comp_bond.";

    const std::vector<std::string> atom_site_attributes = {
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_comp_id",
        "label_seq_id",
        "label_asym_id",
        "label_entity_id",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
    };

    const std::vector<std::string> chem_comp_attributes = {
        "id",
    };

    const std::vector<std::string> chem_comp_bond_attributes = {
        "comp_id",
        "atom_id_1",
        "atom_id_2",
        "value_order",
    };


    for (const auto& molecule : ms.molecules()) {
        std::filesystem::path out_dir{config::chg_out_dir};
        std::string molecule_name = to_lowercase(molecule.name());    
        std::string out_filename = molecule_name + ".fw2.cif";
        std::string out_file{(out_dir / out_filename).string()};
        std::ofstream out_stream{out_file};

        auto document = gemmi::cif::Document{};
        auto& block = document.add_new_block(molecule_name);

        // _atom_site
        auto& atom_site_loop = block.init_loop(atom_site_prefix, atom_site_attributes);
        for (const auto &atom: molecule.atoms()) {
            const std::string group_PDB = atom.hetatm() ? "HETATM" : "ATOM";
            const std::string id = fmt::format("{}", atom.index() + 1);
            const std::string type_symbol = atom.element().symbol();
            const std::string& label_atom_id = id;
            const std::string label_comp_id = atom.residue();
            const std::string label_seq_id = fmt::format("{}", atom.residue_id());
            const std::string label_asym_id = atom.chain_id().empty() ? "." : atom.chain_id();
            const std::string label_entity_id = "1";
            const std::string cartn_x = fmt::format("{:.3f}", atom.pos()[0]);
            const std::string cartn_y = fmt::format("{:.3f}", atom.pos()[1]);
            const std::string cartn_z = fmt::format("{:.3f}", atom.pos()[2]);
            atom_site_loop.add_row({
                group_PDB,
                id,
                type_symbol,
                label_atom_id,
                label_comp_id,
                label_seq_id,
                label_asym_id,
                label_entity_id,
                cartn_x,
                cartn_y,
                cartn_z,
            });
        }

        // _chem_comp
        auto& chem_comp_loop = block.init_loop(chem_comp_prefix, chem_comp_attributes);
        const std::string comp_id = "UNL";
        chem_comp_loop.add_row({
            comp_id,
        });

        // _chem_comp_bond
        auto& chem_comp_bond_loop = block.init_loop(chem_comp_bond_prefix, chem_comp_bond_attributes);
        for (const auto &bond: molecule.bonds()) {
            const std::string atom_id_1 = fmt::format("{}", bond.first().index() + 1);
            const std::string atom_id_2 = fmt::format("{}", bond.second().index() + 1);
            const std::string value_order = convert_bond_order_to_mmcif_value_order_string(bond.order());
            chem_comp_bond_loop.add_row({
                comp_id,
                atom_id_1,
                atom_id_2,
                value_order,
            });
        }

        append_charges_to_block(molecule, charges, block);

        gemmi::cif::write_cif_block_to_stream(out_stream, block);
    }
}

void CIF::save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    auto ext = std::filesystem::path(filename).extension().string();
    
    if (ext == ".cif") {
        generate_mmcif_from_cif_file(ms, charges, config::input_file);
    } else if (ext == ".pdb" or ext == ".ent") {
        generate_mmcif_from_pdb_file(ms, charges, config::input_file);
    } else if (ext == ".mol2" or ext == ".sdf") {
        generate_mmcif_from_atom_and_bond_data(ms, charges);
    }
}
