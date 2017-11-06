set	(SOURCE
	 test_complex_operators.c
	)

set	(LIBS_DEPEND
	 Test_Complex_Containers
	)

set (LIB Test_Complex_Operators)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
