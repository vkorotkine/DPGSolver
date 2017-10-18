add_subdirectory(flux)
add_subdirectory(solution)

set	(SOURCE
	 test_case.c
	)

set	(LIBS_DEPEND
	 Flux
	 Solution
	)

add_library(Test_Case ${SOURCE})
target_link_libraries(Test_Case PRIVATE ${LIBS_DEPEND})
