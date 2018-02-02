set	(SOURCE
	 test_complex_flux.c
	 test_complex_flux_advection.c
	 test_complex_flux_diffusion.c
	 test_complex_flux_euler.c
	 test_complex_flux_navier_stokes.c
	)

set	(LIBS_DEPEND
	 Test_Complex_Containers
	)

set (LIB Test_Complex_Flux)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
