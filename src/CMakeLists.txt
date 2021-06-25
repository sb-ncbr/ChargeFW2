find_package(Boost COMPONENTS program_options system filesystem REQUIRED)
find_package(fmt 5.3 REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(nanoflann 1.3 REQUIRED)
find_package(pybind11)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

configure_file(chargefw2.h.in chargefw2.h)

# Apple Clang does not support -fopenmp; must use llvm clang
# $ brew install llvm libomp
set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY omp)
find_package(OpenMP REQUIRED)

SET(COMMON_LIBS parameters geometry element method config)
foreach (lib ${COMMON_LIBS})
    add_library(${lib} ${lib}.h ${lib}.cpp)
    target_link_libraries(${lib} structures fmt::fmt OpenMP::OpenMP_CXX)
endforeach ()

target_link_libraries(method Eigen3::Eigen)

add_subdirectory(structures)
add_subdirectory(methods)
add_subdirectory(formats)
add_subdirectory(utility)
add_subdirectory(external/lhs)

SET(SOURCES main.cpp periodic_table.cpp periodic_table.h charges.h charges.cpp statistics.h statistics.cpp
        candidates.h candidates.cpp)

add_executable(chargefw2 ${SOURCES})

add_dependencies(chargefw2 ${METHODS})

target_link_libraries(chargefw2 structures ${COMMON_LIBS} formats utility lhs Boost::program_options dl Boost::filesystem fmt::fmt OpenMP::OpenMP_CXX)


set_target_properties(chargefw2 PROPERTIES LINK_FLAGS_RELEASE -s)

install(TARGETS chargefw2 DESTINATION bin)

# Install files required for independent building of new methods
add_library(libchargefw2 STATIC element.h element.cpp method.cpp method.h parameters.h parameters.cpp
        structures/molecule.h structures/molecule.cpp geometry.h geometry.cpp periodic_table.h periodic_table.cpp)

set_target_properties(libchargefw2 PROPERTIES OUTPUT_NAME chargefw2 INTERPROCEDURAL_OPTIMIZATION FALSE)

install(TARGETS libchargefw2 ARCHIVE DESTINATION lib)
install(FILES structures/atom.h structures/bond.h structures/molecule.h DESTINATION include/structures)
install(FILES element.h method.h parameters.h geometry.h periodic_table.h DESTINATION include)

pybind11_add_module(chargefw2_gp python_gp.cpp ${SOURCES})
target_link_libraries(chargefw2_gp PRIVATE libchargefw2 formats ${COMMON_LIBS} utility Boost::program_options)
install(TARGETS chargefw2_gp DESTINATION lib)

pybind11_add_module(chargefw2_python python.cpp ${SOURCES})
target_link_libraries(chargefw2_python PRIVATE libchargefw2 formats ${COMMON_LIBS} utility Boost::program_options)
install(TARGETS chargefw2_python DESTINATION lib)