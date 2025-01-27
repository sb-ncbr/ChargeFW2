set(DEPS method parameters element geometry structures exceptions)

SET(METHODS dummy formal peoe eem sfkeem mgc mpeoe gdac veem kcm denr delre
            tsef charge2 qeq smpqeq eqeq eqeqc abeem sqe sqeq0 sqeqp)
SET(METHODS ${METHODS} PARENT_SCOPE)

foreach(method ${METHODS})
    add_library(${method} SHARED ${method}.cpp ${method}.h)
    target_link_libraries(${method} ${DEPS})
    set_target_properties(${method} PROPERTIES LINK_FLAGS_RELEASE -s SUFFIX ".so")
endforeach(method)

set(NEED_EIGEN eem sfkeem mgc kcm denr delre tsef qeq eqeq eqeqc abeem smpqeq sqe sqeq0 sqeqp)
foreach(method ${NEED_EIGEN})
    target_link_libraries(${method} Eigen3::Eigen)
endforeach(method)

install(TARGETS ${METHODS} LIBRARY DESTINATION lib)
