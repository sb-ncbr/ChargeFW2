cmake_minimum_required(VERSION 3.9)
project(ChargeFW2)

set(CMAKE_CXX_STANDARD 17)


set(MKLROOT /opt/intel/mkl)
set(MKL_CXX_FLAGS  "-DMKL_ILP64 -m64 -I${MKLROOT}/include")
set(MKL_LINKER_FLAGS "-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC -march=native -pedantic ${MKL_CXX_FLAGS}")


file(COPY data/pte.csv DESTINATION data)
add_subdirectory(src)
