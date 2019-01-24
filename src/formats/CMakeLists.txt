add_library(formats sdf.cpp sdf.h reader.h mol2.h mol2.cpp)
target_link_libraries(formats structures fmt::fmt-header-only)