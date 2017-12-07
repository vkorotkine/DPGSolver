set	(SOURCE
	 test_complex_computational_elements.c
	)

set	(LIBS_DEPEND
	)

set(LIB Test_Complex_Computational_Elements)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
