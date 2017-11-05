add_subdirectory(dg)
add_subdirectory(numerical_flux)

set	(SOURCE
	 solve.c
	 solve_explicit.c
	 solve_implicit.c
	)

set	(LIBS_DEPEND
	 ${PETSC_LIBRARIES}
	 Computational_Elements
	 Numerical_Flux
	 Solution
	 Solver_DG
	)

set (LIB Solvers)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
