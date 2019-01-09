find_package(Boost COMPONENTS program_options system filesystem REQUIRED)
find_package(fmt REQUIRED)

include_directories(${PROJECT_BINARY_DIR}/src)

configure_file(config.h.in config.h)

SET(COMMON_LIBS parameters geometry element method)
foreach(lib ${COMMON_LIBS})
    add_library(${lib} ${lib}.h ${lib}.cpp)
    target_link_libraries(${lib} fmt::fmt-header-only)
endforeach()

add_subdirectory(structures)
add_subdirectory(methods)
add_subdirectory(formats)
add_subdirectory(utility)
add_subdirectory(external/lhs)

add_executable(chargefw2 main.cpp periodic_table.cpp periodic_table.h classifier.cpp classifier.h
        charges.cpp charges.h statistics.cpp statistics.h parameterization.h parameterization.cpp utility/utility.cpp)

add_dependencies(chargefw2 ${METHODS})

target_link_libraries(chargefw2 parameters geometry element method formats structures utility lhs
        Boost::program_options Boost::system Boost::filesystem dl stdc++fs nlopt fmt::fmt-header-only)

install(TARGETS chargefw2 DESTINATION bin)
