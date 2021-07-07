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


void CIF::append_fw2_config(gemmi::cif::Block &block) {

    std::string config_prefix = "_chargeFW2_config.";

    std::vector<std::string> config_tags{
        "method", "parameter_set", 
        "read_hetam", "ignore_waters", 
        "permissive_types", "ref_charge_file"
    };

    std::string method = "?";
    std::string param_file = "?";
    std::string read_hetatm = "False";
    std::string ignore_water = "False";
    std::string permissive_types = "False";
    std::string ref_chg_file = "?";

    if (!config::method_name.empty()) method = config::method_name;
    if (!config::par_file.empty()) param_file = fs::path(config::par_file).stem().string();
    if (config::read_hetatm) read_hetatm = "True";
    if (config::ignore_water) ignore_water = "True";
    if (config::permissive_types) permissive_types = "True";
    if (!config::ref_chg_file.empty()) ref_chg_file = fs::path(config::ref_chg_file).filename().string();

    std::vector<std::string> config_data{
        method, param_file,
        read_hetatm, ignore_water,
        permissive_types, ref_chg_file
    };

    for (unsigned i = 0; i != config_tags.size(); ++i) {
            block.set_pair(config_prefix + config_tags[i], config_data[i]);
    }
}


void CIF::replace_fw2_columns(gemmi::cif::Table &table, 
                              std::vector<std::string> &p_charge, 
                              std::vector<std::string> &vdw_radii,
                              const std::vector<std::string> &fw2_tags) {

    std::vector<std::vector<std::string>> fw2_columns = {p_charge, vdw_radii};

    assert(fw2_columns.size() == fw2_tags.size());

    for (unsigned i = 0; i != fw2_tags.size(); ++i){
        auto column = table.bloc.find_loop(fw2_tags[i]);
        std::copy(fw2_columns[i].begin(), fw2_columns[i].end(), column.begin());
    }
}


void CIF::append_fw2_columns(gemmi::cif::Table &table,
                             std::vector<std::string> &p_charge, 
                             std::vector<std::string> &vdw_radii,
                             const std::vector<std::string> &fw2_tags) {

    auto &loop = *table.get_loop();

    unsigned long orig_tag_size = loop.tags.size();
    unsigned long new_tag_size = orig_tag_size + 2;

    // Creates a new table full of empty strings with the correct number of dimensions
    // Outside vector size is the # of columns, inside vector size is the # of rows.
    std::vector<std::vector<std::string>> new_columns(new_tag_size, {loop.length(), {"Empty"}});

    // Copies data from original columns to their respecitve column in the new table filled with empty strings.
    // Leaving only the new appended columns as empty strings
    for (unsigned i = 0; i != orig_tag_size; ++i) {
        auto column = table.bloc.find_loop(loop.tags[i]);
        std::copy(column.begin(), column.end(), new_columns[i].begin());
    }

    new_columns[new_tag_size - 2] = std::move(p_charge);
    new_columns[new_tag_size - 1] = std::move(vdw_radii);

    for (const auto &tag : fw2_tags)
        loop.tags.push_back(tag);

    loop.set_all_values(new_columns);
}                                


void CIF::write_cif_block(std::ostream &out,
                          gemmi::cif::Table &table, 
                          std::vector<std::string> &p_charge, 
                          std::vector<std::string> &vdw_radii) {

    append_fw2_config(table.bloc);

    std::vector<std::string> fw2_tags{
        "_atom_site.fw2_charge",
        "_atom_site.fw2_vdw_radius"};

    auto &loop = *table.get_loop();

    if (loop.has_tag(fw2_tags[0]) && loop.has_tag(fw2_tags[1])){
        replace_fw2_columns(table, p_charge, vdw_radii, fw2_tags);
    } else {
        append_fw2_columns(table, p_charge, vdw_radii, fw2_tags);
    }

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
