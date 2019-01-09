add_library(formats sdf.cpp sdf.h reader.h)
target_link_libraries(formats fmt::fmt-header-only)