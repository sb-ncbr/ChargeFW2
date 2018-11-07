find_package(Qt5 COMPONENTS Core REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Qt5Core_INCLUDE_DIRS})

add_subdirectory(structures)
add_subdirectory(methods)
add_subdirectory(formats)

set(CMAKE_AUTOMOC ON)

add_executable(chargefw2 main.cpp Element.h PeriodicTable.cpp PeriodicTable.h Method.h Classifier.cpp Classifier.h Parameters.cpp Parameters.h Method.cpp)
target_link_libraries(chargefw2 formats methods structures Qt5::Core Eigen3::Eigen)