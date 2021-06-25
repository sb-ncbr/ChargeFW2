cmake_minimum_required(VERSION 3.9)
project(ChargeFW2)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    set(WARNING_FLAGS "-Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded\
        -Wno-exit-time-destructors -Wno-global-constructors -Wno-reserved-id-macro -Wno-sign-conversion\
        -Wno-implicit-int-float-conversion")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(WARNING_FLAGS "-pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wdisabled-optimization\
        -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept\
        -Woverloaded-virtual -Wsign-promo -Wstrict-null-sentinel -Wundef -Werror -Wno-unused -Wno-format-y2k")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_FLAGS} -fPIC -march=native -fopenmp -fsized-deallocation")

# Remove executable permission on shared libraries
set(CMAKE_INSTALL_SO_NO_EXE 1)

# Use LTO for release
if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Make debug builds much faster
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og")

add_subdirectory(src)
add_subdirectory(data)
