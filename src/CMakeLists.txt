add_subdirectory(structures)
add_subdirectory(formats)

add_executable(chargefw2 main.cpp Element.h)
target_link_libraries(chargefw2 formats structures)