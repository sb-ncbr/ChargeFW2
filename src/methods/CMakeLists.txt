set(DEPS method parameters element geometry structures ee_method)

SET(METHODS dummy formal peoe eem sfkeem mgc mpeoe gdac veem kcm denr delre
            tsef charge2 qeq smpqeq eqeq eqeqc abeem sqe sqeq0 sqeqp)

set(METHOD_SOURCES "")
foreach(m ${METHODS})
    list(APPEND METHOD_SOURCES ${m}.cpp ${m}.h)
endforeach()

add_library(methods OBJECT ${METHOD_SOURCES})
target_link_libraries(methods PUBLIC Eigen3::Eigen ${DEPS})
target_include_directories(methods PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
