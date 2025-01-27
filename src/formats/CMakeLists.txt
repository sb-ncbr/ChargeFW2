add_library(formats sdf.cpp sdf.h reader.h mol2.h mol2.cpp pdb.h pdb.cpp mmcif.h mmcif.cpp bonds.cpp bonds.h reader.cpp
            writer.cpp writer.h txt.cpp txt.h pqr.cpp pqr.h common.h common.cpp cif.h cif.cpp)
target_link_libraries(formats structures exceptions fmt::fmt)
