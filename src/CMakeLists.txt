find_package(Boost COMPONENTS program_options system filesystem REQUIRED)

include_directories(${PROJECT_BINARY_DIR}/src)


add_subdirectory(structures)
add_subdirectory(methods)
add_subdirectory(formats)

configure_file(config.h.in config.h)

add_executable(chargefw2 main.cpp element.h periodic_table.cpp periodic_table.h method.h classifier.cpp classifier.h
               parameters.cpp parameters.h method.cpp charges.cpp charges.h element.cpp)
add_dependencies(chargefw2 dummy eem formal peoe sfkeem mgc mpeoe veem kcm denr delre tsef)

target_link_libraries(chargefw2 formats structures Boost::program_options Boost::system Boost::filesystem dl stdc++fs)

install(TARGETS chargefw2 DESTINATION bin)
