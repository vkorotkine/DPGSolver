set	(SOURCE
	 test_complex_compute_all_rhs_dpg.c
	 test_complex_solve_dpg.c
	 test_complex_face_solver_dpg.c
	 test_complex_volume_solver_dpg.c
	)

set	(LIBS_DEPEND
	 Test_Support_Containers
	 Test_Support_Computational_Elements
	 Test_Complex_Solvers
	)

set (LIB Test_Complex_Solver_DPG)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
