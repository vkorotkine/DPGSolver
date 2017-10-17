add_subdirectory(dg)

set	(SOURCE
	 solve.c
	 solve_explicit.c
	 solve_implicit.c
	)

set	(LIBS_DEPEND
	 Computational_Elements
	 Solution
	 Solver_DG
	)

set (LIB Solvers)
add_library(${LIB} ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
