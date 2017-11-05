add_subdirectory(containers)
add_subdirectory(simulation)

set	(SOURCE
	)

set	(LIBS_DEPEND
	 Test_Complex_Containers
	 Test_Complex_Flux
	 Test_Complex_Solver_DG
	)

set (LIB Test_Complex)
add_library(${LIB} ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
