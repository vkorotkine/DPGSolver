set	(SOURCE
	 compute_face_rlhs_dg.c
	 compute_grad_coef_dg.c
	 compute_volume_rlhs_dg.c
	 element_solver_dg.c
	 face_solver_dg.c
	 volume_solver_dg.c
	)

set	(LIBS_DEPEND
	 Containers
	 Element
	 Flux
	)

set (LIB Solver_DG)
add_library(${LIB} ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
