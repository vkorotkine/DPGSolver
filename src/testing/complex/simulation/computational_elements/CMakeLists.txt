set	(SOURCE
	 test_complex_computational_elements.c
	)

set	(LIBS_DEPEND
	 Test_Complex_Solvers
	)

set(LIB Test_Complex_Computational_Elements)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
