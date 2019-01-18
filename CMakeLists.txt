cmake_minimum_required(VERSION 3.9)
project(ChargeFW2)

set(CMAKE_CXX_STANDARD 17)

# Setup MKL
set(MKLROOT /opt/intel/mkl)
set(MKL_CXX_FLAGS  "-m64 -isystem ${MKLROOT}/include")
set(MKL_LINKER_FLAGS "-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl")


set(WARNING_FLAGS "-pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization\
                   -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept\
                   -Wold-style-cast -Woverloaded-virtual -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel\
                   -Wstrict-overflow=5 -Wundef -Werror -Wno-unused")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_FLAGS} -fPIC -march=native")


# Remove executable permission on shared libraries
set(CMAKE_INSTALL_SO_NO_EXE 1)

# Strip binaries in release
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s")


# Make debug builds much faster
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og")

add_subdirectory(src)
add_subdirectory(data)