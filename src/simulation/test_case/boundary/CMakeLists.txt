set	(SOURCE
	 boundary.c
	 boundary_advection.c
	 boundary_diffusion.c
	 boundary_euler.c
	 boundary_navier_stokes.c
	)

set	(LIBS_DEPEND
	 Solvers
	)

set (LIB Boundary)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
