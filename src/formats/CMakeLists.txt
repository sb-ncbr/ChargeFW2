add_library(formats sdf.cpp sdf.h reader.h)
target_link_libraries(formats structures fmt::fmt-header-only)