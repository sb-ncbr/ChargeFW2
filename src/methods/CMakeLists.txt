add_library(formal SHARED Formal.cpp Formal.h ../Method.h ../Method.cpp ../Parameters.cpp ../Parameters.h)
add_library(dummy SHARED Dummy.cpp Dummy.h ../Method.h ../Method.cpp ../Parameters.cpp ../Parameters.h)
add_library(eem SHARED EEM.cpp EEM.h ../Method.h ../Method.cpp ../Parameters.cpp ../Parameters.h)

set_target_properties(eem PROPERTIES COMPILE_FLAGS ${MKL_CXX_FLAGS})
target_link_libraries(eem ${MKL_LINKER_FLAGS})

install(TARGETS dummy eem formal LIBRARY DESTINATION lib)