cmake_minimum_required(VERSION 3.9)
project(ChargeFW2)

set(CMAKE_CXX_STANDARD 17)

# Setup MKL
set(MKLROOT /opt/intel/mkl)
set(MKL_CXX_FLAGS  "-m64 -I${MKLROOT}/include")
set(MKL_LINKER_FLAGS "-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -fPIC -march=native")

# Strip binaries in release
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s")

add_subdirectory(src)
add_subdirectory(data)