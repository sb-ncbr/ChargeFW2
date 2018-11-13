find_package(Boost COMPONENTS program_options system filesystem REQUIRED)

add_subdirectory(structures)
add_subdirectory(methods)
add_subdirectory(formats)

add_executable(chargefw2 main.cpp Element.h PeriodicTable.cpp PeriodicTable.h Method.h Classifier.cpp Classifier.h Parameters.cpp Parameters.h Method.cpp Charges.cpp Charges.h)
add_dependencies(chargefw2 dummy eem formal peoe)

target_link_libraries(chargefw2 formats structures Boost::program_options Boost::system Boost::filesystem dl)

install(TARGETS chargefw2 DESTINATION bin)
