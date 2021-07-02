//
// Created by danny305 on May 3rd, 2021.
//

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <filesystem>

#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>

#define GEMMI_WRITE_IMPLEMENTATION
#include <gemmi/to_cif.hpp>
#undef GEMMI_WRITE_IMPLEMENTATION

#include "chargefw2.h"
#include "cif.h"
#include "../structures/molecule_set.h"
#include "../charges.h"
#include "../config.h"


namespace fs = std::filesystem;

static const std::vector<std::string> atom_site_columns{
        "group_PDB",
        "auth_asym_id",
        "auth_seq_id",
        "pdbx_PDB_ins_code",
        "auth_comp_id",
        "auth_atom_id",
        "label_alt_id",
        "type_symbol",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "pdbx_PDB_model_num",
    };


class MCRA {
    const int _model;
    const std::string _chain;
    const std::string _res_num;
    const std::string _residue;
    const std::string _atom;

public:
    MCRA(const int model,
         const std::string &chain,
         const std::string &res_num,
         const std::string &residue,
         const std::string &atom)
        : _model(model), _chain(chain), _res_num(res_num), _residue(residue), _atom(atom)
        {}

    int find_row(gemmi::cif::Table &table, int start_idx = 0) const {
        int idx = start_idx; 
        int total_rows = static_cast<int>(table.length());

        assert(idx < total_rows);

        while (idx < total_rows) {
            if (this->is_row(table[idx])) {
                return idx;
            }
            ++idx;
        }

        fmt::print(stderr, "In MCRA::find_row(), unable to find Atom(1, {}, {}, {}, {})\n"
                    "Looping through entire table to double check...",
                    this->_chain, this->_res_num, this->_residue, this->_atom); 

        idx = 0;
        for (const auto site : table) {
            if (this->is_row(site)) {
                return idx;
            }
            ++idx;
        }
        return -1;
    }

    friend std::ostream &operator<<(std::ostream &os, const MCRA &mcra) {
        os << "Atom(" << mcra._model << " " << mcra._chain << " " << mcra._res_num << " [" << mcra._residue << "] " << mcra._atom << ")";
        return os;
    }

private:
    bool is_row(const gemmi::cif::Table::Row &row) const {
        int model = static_cast<int>(std::stoul(row[11]));
        const std::string &chain = row[1];
        const std::string &residue = row[4];
        const std::string &res_num = row[2];
        const std::string atom = row[5][0] != '"' ? row[5] : row[5].substr(1, row[5].size() - 2);

        if (_model == model) {
            if (_chain == chain) {
                if (_res_num == res_num && _residue == residue) {
                    if (_atom == atom) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
};


void CIF::write_cif_block(std::ostream &out,
                          gemmi::cif::Table &table, 
                          std::vector<std::string> &p_charge, 
                          std::vector<std::string> &vdw_radii) {

    auto &loop = *table.get_loop();

    unsigned long orig_tag_size = loop.tags.size();
    unsigned long new_tag_size = orig_tag_size + 2;

    // Creates a new table full of empty strings with the correct number of dimensions
    // Outside vector size is the # of columns, inside vector size is the # of rows.
    std::vector<std::vector<std::string>> new_columns(new_tag_size, {loop.length(), {"Empty"}});

    // Copies data from original columns to their respecitve column in the new table filled with empty strings.
    // Leaving only the new appended columns as empty strings
    for (unsigned int i = 0; i != orig_tag_size; ++i) {
        auto column = table.bloc.find_loop(loop.tags[i]);
        std::copy(column.begin(), column.end(), new_columns[i].begin());
    }

    new_columns[new_tag_size - 2] = std::move(p_charge);
    new_columns[new_tag_size - 1] = std::move(vdw_radii);

    std::vector<std::string> new_tags{
        "_atom_site.fw2_charge",
        "_atom_site.fw2_vdw_radius"};
    for (const auto &tag : new_tags)
        loop.tags.push_back(tag);

    loop.set_all_values(new_columns);

    gemmi::cif::write_cif_block_to_stream(out, table.bloc);
}


void CIF::save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {

    fs::path out_dir{config::chg_out_dir};
    std::string out_filename = fs::path(filename).filename().replace_extension(".fw2.cif").string();
    std::string out_file{(out_dir / out_filename).string()};
    std::ofstream out_stream{out_file};

    auto doc = gemmi::cif::read_file(filename);
    auto& block = doc.sole_block();
    auto table = block.find("_atom_site.", atom_site_columns);

    if (!table.ok()) {
        fmt::print(stderr, "_atom_site category is empty for file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    auto p_charge  = std::vector<std::string>{table.length(), "?"};
    auto vdw_radii = std::vector<std::string>{table.length(), "?"};

    const auto &molecule = ms.molecules()[0];


    // ChargeFW2 is hardcoded to only read first model.
    const int model = 1;
    int row_num = 0;

    try {
        auto chg = charges[molecule.name()];
        for (size_t i = 0; i < molecule.atoms().size(); i++) {
            const auto &atom = molecule.atoms()[i];

            MCRA mcra{
                model, 
                atom.chain_id(), 
                std::to_string(atom.residue_id()), 
                atom.residue(), 
                atom.name()
            };

            row_num = mcra.find_row(table, row_num);

            if (row_num == -1){
                 fmt::print(stderr, "Failed to find Atom(1, {}, {}, {}, {})\n", 
                    atom.chain_id(), atom.residue_id(), atom.residue(), atom.name());
                exit(EXIT_FILE_ERROR);
            }
            p_charge[row_num]  = std::to_string(chg[i]);
            vdw_radii[row_num] = std::to_string(atom.element().vdw_radius());
        }

        write_cif_block(out_stream, table, p_charge, vdw_radii);
    }
    catch (std::out_of_range &) {
        /* Do nothing */
    }
}
