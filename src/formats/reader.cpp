//
// Created by krab1k on 29.1.19.
//

#include <filesystem>

#include "reader.h"
#include "sdf.h"
#include "pdb.h"
#include "mmcif.h"
#include "mol2.h"
#include "../utility/strings.h"


Reader::~Reader() = default;

MoleculeSet load_molecule_set(const std::string &filename) {
    auto ext = std::filesystem::path(filename).extension().string();

    std::unique_ptr<Reader> reader;
    ext = to_lowercase(ext);
    if (ext == ".sdf") {
        reader = std::make_unique<SDF>();
    } else if (ext == ".mol2") {
        reader = std::make_unique<Mol2>();
    } else if (ext == ".pdb" or ext == ".ent") {
        reader = std::make_unique<PDB>();
    } else if (ext == ".cif") {
        reader = std::make_unique<mmCIF>();
    } else {
        fmt::print(stderr, "Filetype {} not supported\n", ext);
        exit(EXIT_FILE_ERROR);
    }

    return reader->read_file(filename);
}
