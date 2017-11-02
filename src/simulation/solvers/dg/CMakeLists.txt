set	(SOURCE
	 complex/volume_solver_dg_complex.c
	 compute_face_rlhs_dg.c
	 compute_grad_coef_dg.c
	 compute_source_rlhs_dg.c
	 compute_volume_rlhs_dg.c
	 element_solver_dg.c
	 face_solver_dg.c
	 solve_dg.c
	 volume_solver_dg.c
	)

set	(LIBS_DEPEND
	 ${PETSC_LIBRARIES}
	 Complex_Containers
	 Containers
	 Element
	 Flux
	 Numerical_Flux
	)

set (LIB Solver_DG)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
