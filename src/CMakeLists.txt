find_package(Boost COMPONENTS program_options system filesystem REQUIRED)

add_subdirectory(structures)
add_subdirectory(methods)
add_subdirectory(formats)

configure_file(config.h.in config.h)

add_executable(chargefw2 main.cpp element.h periodic_table.cpp periodic_table.h method.h classifier.cpp classifier.h parameters.cpp parameters.h method.cpp charges.cpp charges.h config.h.in)
add_dependencies(chargefw2 dummy eem formal peoe)

target_include_directories(chargefw2 PRIVATE ${PROJECT_BINARY_DIR}/src)
target_link_libraries(chargefw2 formats structures Boost::program_options Boost::system Boost::filesystem dl)

install(TARGETS chargefw2 DESTINATION bin)
