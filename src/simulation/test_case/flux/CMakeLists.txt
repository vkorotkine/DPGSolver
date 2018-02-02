set	(SOURCE
	 flux.c
	 flux_advection.c
	 flux_diffusion.c
	 flux_euler.c
	 flux_navier_stokes.c
	)

set	(LIBS_DEPEND
	 Containers
	)

set (LIB Flux)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
