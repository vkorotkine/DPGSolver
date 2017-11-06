add_subdirectory(boundary)
add_subdirectory(flux)

set	(SOURCE
	 test_complex_test_case.c
	)

set	(LIBS_DEPEND
	)

set (LIB Test_Complex_Test_Case)
add_library(${LIB} ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
