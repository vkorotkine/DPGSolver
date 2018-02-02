set	(SOURCE
	 compute_error.c
	 compute_error_advection.c
	 compute_error_diffusion.c
	 compute_error_euler.c
	 compute_error_navier_stokes.c
	)

set	(LIBS_DEPEND
	 ${MPI_C_LIBRARIES}
	 Computational_Elements
	 Containers
	 Element
	 Solvers
	)

add_library(Compute_Error STATIC ${SOURCE})
target_link_libraries(Compute_Error ${LIBS_DEPEND})
