cmake_minimum_required(VERSION 3.9)
project(ChargeFW2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC -march=native -pedantic")

file(COPY data/pte.csv DESTINATION data)
add_subdirectory(src)
