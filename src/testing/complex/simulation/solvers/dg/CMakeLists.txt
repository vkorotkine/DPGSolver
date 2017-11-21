set	(SOURCE
	 test_complex_compute_face_rhs_dg.c
	 test_complex_compute_volume_rhs_dg.c
	 test_complex_solve_dg.c
	)

set	(LIBS_DEPEND
	 Test_Complex_Boundary
	 Test_Complex_Containers
	 Test_Complex_Flux
	 Test_Complex_Numerical_Flux
	 Test_Complex_Operators
	 Test_Complex_Test_Case
	 Test_Support_Computational_Elements
	 Test_Support_General
	)

set (LIB Test_Complex_Solver_DG)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
