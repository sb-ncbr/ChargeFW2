add_library(dummy SHARED Dummy.cpp Dummy.h ../Method.h ../Method.cpp ../Parameters.cpp ../Parameters.h)
add_library(eem SHARED EEM.cpp EEM.h ../Method.h ../Method.cpp ../Parameters.cpp ../Parameters.h)

install(TARGETS dummy eem LIBRARY DESTINATION lib)