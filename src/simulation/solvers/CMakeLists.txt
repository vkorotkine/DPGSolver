add_subdirectory(dg)
add_subdirectory(dpg)
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
	)

set	(LIBS_DEPEND
	 ${PETSC_LIBRARIES}
	 Complex_Containers
#	 Computational_Elements
#	 Numerical_Flux
#	 Solution
	 Solver_DG
	 Solver_DPG
	)

set (LIB Solvers)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
