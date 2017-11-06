set	(SOURCE
	 test_complex_flux.c
	 test_complex_flux_advection.c
	)

set	(LIBS_DEPEND
	 Test_Complex_Containers
	)

set (LIB Test_Complex_Flux)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
