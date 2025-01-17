find_package(Boost COMPONENTS program_options REQUIRED)
find_package(fmt 5.3 REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(nanoflann 1.3 REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(gemmi 0.6.2 CONFIG REQUIRED)

option(PYTHON_MODULE "Build Python bindings" ON)

if(PYTHON_MODULE)
    find_package(pybind11)
endif ()

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

configure_file(chargefw2.h.in chargefw2.h)

SET(COMMON_LIBS parameters geometry element method config)
foreach (lib ${COMMON_LIBS})
    add_library(${lib} ${lib}.h ${lib}.cpp)
    target_link_libraries(${lib} structures exceptions fmt::fmt)
endforeach ()

target_link_libraries(method Eigen3::Eigen)

add_subdirectory(structures)
add_subdirectory(methods)
add_subdirectory(formats)
add_subdirectory(exceptions)
add_subdirectory(utility)

SET(SOURCES periodic_table.cpp periodic_table.h charges.h charges.cpp candidates.h candidates.cpp)

add_library(common ${SOURCES})

add_executable(chargefw2 main.cpp options.cpp options.h)
add_dependencies(chargefw2 ${METHODS})

target_link_libraries(chargefw2 structures ${COMMON_LIBS} formats common exceptions utility Boost::program_options dl fmt::fmt gemmi::gemmi_cpp)

set_target_properties(chargefw2 PROPERTIES LINK_FLAGS_RELEASE -s)

install(TARGETS chargefw2 DESTINATION bin)

if(PYTHON_MODULE)
    pybind11_add_module(chargefw2_python python.cpp)
    target_link_libraries(chargefw2_python PRIVATE formats common ${COMMON_LIBS} utility gemmi::gemmi_cpp)
    set_target_properties(chargefw2_python PROPERTIES OUTPUT_NAME chargefw2)
    install(TARGETS chargefw2_python DESTINATION lib)
endif ()
