set	(SOURCE
	 test_complex_boundary.c
	 test_complex_boundary_advection.c
	 test_complex_boundary_euler.c
	)

set	(LIBS_DEPEND
	 Complex_Containers
	 Test_Complex_Containers
	 Test_Complex_Operators
	 Test_Complex_Solver_DG
	)

set (LIB Test_Complex_Boundary)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
