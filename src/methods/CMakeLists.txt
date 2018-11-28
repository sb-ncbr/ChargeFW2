add_library(dummy SHARED dummy.cpp dummy.h ../method.h ../method.cpp ../parameters.cpp ../parameters.h)
add_library(formal SHARED formal.cpp formal.h ../method.h ../method.cpp ../parameters.cpp ../parameters.h)
add_library(peoe SHARED peoe.cpp peoe.h ../method.h ../method.cpp ../parameters.cpp ../parameters.h)
add_library(eem SHARED eem.cpp eem.h ../method.h ../method.cpp ../parameters.cpp ../parameters.h)

set_target_properties(eem PROPERTIES COMPILE_FLAGS ${MKL_CXX_FLAGS})
target_link_libraries(eem ${MKL_LINKER_FLAGS})

install(TARGETS dummy formal peoe eem LIBRARY DESTINATION lib)
