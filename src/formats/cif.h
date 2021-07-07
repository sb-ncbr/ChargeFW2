//
// Created by danny305 on May 3rd 2021.
//

#pragma once

#include <string>
#include <vector>
#include <gemmi/cif.hpp>

#include "writer.h"
#include "../structures/molecule_set.h"
#include "../charges.h"


class CIF : public Writer {
public:
    void save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) override;
private:
    static void append_fw2_config(gemmi::cif::Block &block, const std::string &method, const std::string &parameters);
    static void replace_fw2_columns(gemmi::cif::Table &table, const std::vector<std::string> &p_charge,
                                    const std::vector<std::string> &vdw_radii,
                                    const std::vector<std::string> &fw2_tags = {"_atom_site.fw2_charge","_atom_site.fw2_vdw_radius"});
    static void append_fw2_columns(gemmi::cif::Table &table, const std::vector<std::string> &p_charge,
                                   const std::vector<std::string> &vdw_radii,
                                   const std::vector<std::string> &fw2_tags = {"_atom_site.fw2_charge","_atom_site.fw2_vdw_radius"});
    static void write_cif_block(std::ostream &out, gemmi::cif::Table &table,
                                const std::vector<std::string> &p_charge, const std::vector<std::string> &vdw_radii,
                                const std::string &method, const std::string &parameters);
};
