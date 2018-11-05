find_package(Qt5 COMPONENTS Core REQUIRED)
include_directories(${Qt5Core_INCLUDE_DIRS})

add_subdirectory(structures)
add_subdirectory(formats)
set(CMAKE_AUTOMOC ON)

add_executable(chargefw2 main.cpp Element.h PeriodicTable.cpp PeriodicTable.h Method.h methods/EEM.cpp methods/EEM.h Classifier.cpp Classifier.h)
target_link_libraries(chargefw2 formats structures Qt5::Core)