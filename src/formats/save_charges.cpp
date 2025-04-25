#include <exception>
#include <filesystem>
#include <fmt/core.h>

#include "../charges.h"
#include "../config.h"
#include "../structures/molecule_set.h"
#include "cif.h"
#include "mol2.h"
#include "pqr.h"
#include "save_charges.h"
#include "txt.h"

void save_charges(MoleculeSet &ms, Charges &charges,
                  const std::string &filename) {
  std::filesystem::path out_dir(config::chg_out_dir);
  auto file_path = std::filesystem::path(filename);
  auto ext = file_path.extension().string();

  auto txt_str = file_path.filename().string() + ".txt";
  TXT().save_charges(ms, charges, out_dir / std::filesystem::path(txt_str));

  CIF().save_charges(ms, charges, filename);

  if (ms.has_proteins()) {
    auto pqr_str = file_path.filename().string() + ".pqr";
    PQR().save_charges(ms, charges, out_dir / std::filesystem::path(pqr_str));
  } else {
    auto mol2_str = file_path.filename().string() + ".mol2";
    Mol2().save_charges(ms, charges,
                        out_dir / std::filesystem::path(mol2_str));
  }
}
