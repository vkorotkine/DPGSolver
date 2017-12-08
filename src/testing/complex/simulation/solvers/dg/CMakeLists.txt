set	(SOURCE
	 test_complex_compute_face_rhs_dg.c
	 test_complex_compute_volume_rhs_dg.c
	 test_complex_solve_dg.c
	 test_complex_face_solver_dg.c
	 test_complex_volume_solver_dg.c
	)

set	(LIBS_DEPEND
	 Test_Support_Containers
	 Test_Support_Computational_Elements
	 Test_Complex_Solvers
	)

set (LIB Test_Complex_Solver_DG)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
