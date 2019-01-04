find_package(Boost COMPONENTS program_options system filesystem REQUIRED)

include_directories(${PROJECT_BINARY_DIR}/src)

configure_file(config.h.in config.h)

add_library(parameters parameters.h parameters.cpp)
add_library(geometry geometry.h geometry.cpp)
add_library(element element.h element.cpp)
add_library(method method.h method.cpp)

add_subdirectory(structures)
add_subdirectory(methods)
add_subdirectory(formats)

add_executable(chargefw2 main.cpp periodic_table.cpp periodic_table.h classifier.cpp classifier.h
               charges.cpp charges.h statistics.cpp statistics.h)

add_dependencies(chargefw2 ${METHODS})

target_link_libraries(chargefw2 parameters geometry element method formats structures
                      Boost::program_options Boost::system Boost::filesystem dl stdc++fs)

install(TARGETS chargefw2 DESTINATION bin)
