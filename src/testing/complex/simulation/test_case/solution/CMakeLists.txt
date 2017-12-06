set	(SOURCE
	 test_complex_solution.c
	 test_complex_solution_euler.c
	)

set	(LIBS_DEPEND
	)

set(LIB Test_Complex_Solution)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
