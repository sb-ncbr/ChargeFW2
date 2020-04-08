//
// Created by krab1k on 28.1.19.
//


#pragma once

#include <gemmi/cif.hpp>

#include "reader.h"


class mmCIF: public Reader {
    static void process_record(const std::string &structure_data, std::unique_ptr<std::vector<Molecule>> &molecules);

    static void read_protein_molecule(gemmi::cif::Block &data, std::unique_ptr<std::vector<Atom>> &atoms);

    static void read_ccd_molecule(gemmi::cif::Block &data, std::unique_ptr<std::vector<Atom>> &atoms, std::unique_ptr<std::vector<Bond>> &bonds);

public:
    mmCIF();

    MoleculeSet read_file(const std::string &filename) override;
};
