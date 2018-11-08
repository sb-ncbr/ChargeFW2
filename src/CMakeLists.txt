find_package(Boost REQUIRED)

add_subdirectory(structures)
add_subdirectory(methods)
add_subdirectory(formats)

set(CMAKE_AUTOMOC ON)

add_executable(chargefw2 main.cpp Element.h PeriodicTable.cpp PeriodicTable.h Method.h Classifier.cpp Classifier.h Parameters.cpp Parameters.h Method.cpp)
target_link_libraries(chargefw2 formats methods structures ${MKL_LINKER_FLAGS})