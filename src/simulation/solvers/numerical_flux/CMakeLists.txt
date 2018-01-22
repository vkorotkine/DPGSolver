set	(SOURCE
	 numerical_flux.c
	 numerical_flux_advection.c
	 numerical_flux_diffusion.c
	 numerical_flux_euler.c
	)

set	(LIBS_DEPEND
	 Solvers
	)

set (LIB Numerical_Flux)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
