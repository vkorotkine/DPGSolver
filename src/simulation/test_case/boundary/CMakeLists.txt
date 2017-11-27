set	(SOURCE
	 boundary.c
	 boundary_advection.c
	 boundary_euler.c
	)

set	(LIBS_DEPEND
	 Solver_DG
	)

set (LIB Boundary)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
