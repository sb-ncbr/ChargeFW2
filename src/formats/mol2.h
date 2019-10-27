//
// Created by krab1k on 24.1.19.
//

#pragma once

#include "reader.h"
#include "writer.h"
#include "../charges.h"


class Mol2 : public Reader, public Writer {
    static void read_until_end_of_record(std::ifstream &file);

    static void read_record(std::ifstream &file, std::string &line, std::unique_ptr<std::vector<Atom>> &atoms,
                     std::unique_ptr<std::vector<Bond>> &bonds);

public:
    MoleculeSet read_file(const std::string &filename) override;

    void save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) override;
};
