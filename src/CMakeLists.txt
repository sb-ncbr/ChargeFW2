find_package(Boost COMPONENTS program_options system filesystem REQUIRED)
find_package(fmt REQUIRED)
find_package(nanoflann REQUIRED)

include_directories(${PROJECT_BINARY_DIR}/src ${EIGEN3_INCLUDE_DIR})

configure_file(chargefw2.h.in chargefw2.h)

SET(COMMON_LIBS parameters geometry element method config)
foreach (lib ${COMMON_LIBS})
    add_library(${lib} ${lib}.h ${lib}.cpp)
    target_link_libraries(${lib} structures fmt::fmt)
endforeach ()

target_link_libraries(method Eigen3::Eigen)

add_subdirectory(structures)
add_subdirectory(methods)
add_subdirectory(formats)
add_subdirectory(utility)
add_subdirectory(external/lhs)

SET(SOURCES main.cpp periodic_table.cpp periodic_table.h charges.h charges.cpp statistics.h statistics.cpp
        parameterization.h parameterization.cpp candidates.h candidates.cpp)

add_executable(chargefw2 ${SOURCES})

add_dependencies(chargefw2 ${METHODS})

target_link_libraries(chargefw2 structures ${COMMON_LIBS} formats utility lhs
        Boost::program_options Boost::system Boost::filesystem dl stdc++fs nlopt_cxx fmt::fmt)

install(TARGETS chargefw2 DESTINATION bin)

# Install files required for independent building of new methods
add_library(libchargefw2 STATIC element.h element.cpp method.cpp method.h parameters.h parameters.cpp
        structures/molecule.h structures/molecule.cpp geometry.h geometry.cpp periodic_table.h periodic_table.cpp)

set_target_properties(libchargefw2 PROPERTIES OUTPUT_NAME chargefw2)

install(TARGETS libchargefw2 LIBRARY DESTINATION lib)
install(FILES structures/atom.h structures/bond.h structures/molecule.h DESTINATION include/structures)
install(FILES element.h method.h parameters.h geometry.h periodic_table.h DESTINATION include)
