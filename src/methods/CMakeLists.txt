add_library(dummy SHARED dummy.cpp dummy.h ../method.h ../method.cpp ../parameters.cpp ../parameters.h)
add_library(formal SHARED formal.cpp formal.h ../method.h ../method.cpp ../parameters.cpp ../parameters.h)
add_library(peoe SHARED peoe.cpp peoe.h ../method.h ../method.cpp ../parameters.cpp ../parameters.h)
add_library(mpeoe SHARED mpeoe.cpp mpeoe.h ../method.h ../method.cpp ../parameters.cpp ../parameters.h)
add_library(gdac SHARED gdac.cpp gdac.h ../method.h ../method.cpp ../parameters.cpp ../parameters.h ../element.h)
add_library(eem SHARED eem.cpp eem.h ../method.h ../method.cpp ../parameters.cpp ../parameters.h)
add_library(sfkeem SHARED sfkeem.cpp sfkeem.h ../method.h ../method.cpp ../parameters.cpp ../parameters.h)
add_library(mgc SHARED mgc.cpp mgc.h ../method.h ../method.cpp ../structures/molecule.cpp)
add_library(veem SHARED veem.cpp veem.h ../method.h ../method.cpp ../structures/molecule.cpp ../element.cpp)

set_target_properties(eem PROPERTIES COMPILE_FLAGS ${MKL_CXX_FLAGS})
target_link_libraries(eem ${MKL_LINKER_FLAGS})

set_target_properties(sfkeem PROPERTIES COMPILE_FLAGS ${MKL_CXX_FLAGS})
target_link_libraries(sfkeem ${MKL_LINKER_FLAGS})

set_target_properties(mgc PROPERTIES COMPILE_FLAGS ${MKL_CXX_FLAGS})
target_link_libraries(mgc ${MKL_LINKER_FLAGS})

install(TARGETS dummy formal peoe eem sfkeem mgc mpeoe gdac veem LIBRARY DESTINATION lib)
