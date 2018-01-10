add_subdirectory(adaptation)
add_subdirectory(numerical_flux)

set	(SOURCE
	 compute_face_rlhs.c
	 compute_volume_rlhs.c
	 element_solver.c
	 face_solver.c
	 solve.c
	 solve_explicit.c
	 solve_implicit.c
	 volume_solver.c

	 dg/compute_face_rlhs_dg.c
	 dg/compute_grad_coef_dg.c
	 dg/compute_source_rlhs_dg.c
	 dg/compute_volume_rlhs_dg.c
	 dg/element_solver_dg.c
	 dg/face_solver_dg.c
	 dg/solve_dg.c
	 dg/volume_solver_dg.c

	 dpg/compute_all_rlhs_dpg.c
	 dpg/face_solver_dpg.c
	 dpg/element_solver_dpg.c
	 dpg/solve_dpg.c
	 dpg/volume_solver_dpg.c
	)

set	(LIBS_DEPEND
	 ${PETSC_LIBRARIES}
	 Boundary
	 Post_Processing
	 Test_Complex_Computational_Elements
	 Test_Complex_Test_Case
	)

set (LIB Solvers)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
