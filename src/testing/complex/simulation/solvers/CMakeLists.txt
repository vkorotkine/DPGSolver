add_subdirectory(numerical_flux)

set	(SOURCE
	 test_complex_compute_face_rhs.c
	 test_complex_compute_volume_rhs.c
	 test_complex_face_solver.c
	 test_complex_solve.c
	 test_complex_volume_solver.c

	 dg/test_complex_compute_face_rhs_dg.c
	 dg/test_complex_compute_volume_rhs_dg.c
	 dg/test_complex_solve_dg.c
	 dg/test_complex_face_solver_dg.c
	 dg/test_complex_volume_solver_dg.c

	 dpg/test_complex_compute_all_rhs_dpg.c
	 dpg/test_complex_solve_dpg.c
	 dpg/test_complex_face_solver_dpg.c
	 dpg/test_complex_volume_solver_dpg.c
	)

set	(LIBS_DEPEND
	 Test_Complex_Boundary
	 Test_Complex_Flux
	 Test_Complex_Geometry
	 Test_Support_Containers
	 Test_Support_Computational_Elements
	)

set (LIB Test_Complex_Solvers)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
