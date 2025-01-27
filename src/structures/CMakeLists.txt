add_library(structures atom.cpp atom.h molecule.cpp molecule.h molecule_set.cpp molecule_set.h bond.h bond.cpp)
target_link_libraries(structures exceptions geometry utility fmt::fmt)
