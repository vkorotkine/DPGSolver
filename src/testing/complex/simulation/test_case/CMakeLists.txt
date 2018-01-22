add_subdirectory(boundary)
add_subdirectory(flux)
add_subdirectory(solution)

set	(SOURCE
	 test_complex_test_case.c
	)

set	(LIBS_DEPEND
	 Test_Case
	 Test_Complex_Boundary
	 Test_Complex_Solution
	)

set (LIB Test_Complex_Test_Case)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
