set	(SOURCE
	 complex/face_solver_dpg_complex.c
	 complex/volume_solver_dpg_complex.c
	 compute_all_rlhs_dpg.c
	 face_solver_dpg.c
	 element_solver_dpg.c
	 solve_dpg.c
	 volume_solver_dpg.c
	)

set	(LIBS_DEPEND
#	 ${PETSC_LIBRARIES}
	 Solvers
	)

set (LIB Solver_DPG)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
