set(DEPS method parameters element geometry structures)

SET(METHODS dummy formal peoe eem sfkeem mgc mpeoe gdac veem kcm denr delre tsef charge2 qeq smpqeq eqeq eqeqc abeem)
SET(METHODS ${METHODS} PARENT_SCOPE)

foreach(method ${METHODS})
    add_library(${method} SHARED ${method}.cpp ${method}.h)
    target_link_libraries(${method} ${DEPS})
endforeach(method)

set(NEED_MKL eem sfkeem mgc kcm denr delre tsef qeq eqeq eqeqc abeem)
foreach(method ${NEED_MKL})
    set_target_properties(${method} PROPERTIES COMPILE_FLAGS ${MKL_CXX_FLAGS})
    target_link_libraries(${method} ${MKL_LINKER_FLAGS})
endforeach(method)

install(TARGETS ${METHODS} LIBRARY DESTINATION lib)
