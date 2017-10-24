add_subdirectory(flux)
add_subdirectory(solution)
add_subdirectory(compute_error)

set	(SOURCE
	 test_case.c
	)

set	(LIBS_DEPEND
	 Compute_Error
	 Flux
	 Numerical_Flux
	 Solution
	)

add_library(Test_Case ${SOURCE})
target_link_libraries(Test_Case PRIVATE ${LIBS_DEPEND})
