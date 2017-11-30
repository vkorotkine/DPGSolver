set	(SOURCE
	 test_complex_numerical_flux.c
	 test_complex_numerical_flux_advection.c
	 test_complex_numerical_flux_euler.c
	)

set	(LIBS_DEPEND
	 Numerical_Flux
	 Test_Complex_General
	 Test_Complex_Containers
	)

set (LIB Test_Complex_Numerical_Flux)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
