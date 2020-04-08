//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <string>

#include "reader.h"


class SDF: public Reader {
    static void read_V2000(std::ifstream &file, std::string &line, std::unique_ptr<std::vector<Atom>> &atoms,
                    std::unique_ptr<std::vector<Bond>> &bonds);

    static void read_V3000(std::ifstream &file, std::string &line, std::unique_ptr<std::vector<Atom>> &atoms,
                    std::unique_ptr<std::vector<Bond>> &bonds);

    static void read_until_end_of_record(std::ifstream &file);

public:
    SDF();
    MoleculeSet read_file(const std::string &filename) override;
};
