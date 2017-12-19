set	(SOURCE
	 compute_all_rlhs_dpg.c
	 face_solver_dpg.c
	 element_solver_dpg.c
	 solve_dpg.c
	 volume_solver_dpg.c
	)

set	(LIBS_DEPEND
	 Solvers
	 Test_Complex_Computational_Elements
	 Test_Complex_Test_Case
	)

set (LIB Solver_DPG)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
