set	(SOURCE
	 numerical_flux.c
	 numerical_flux_euler.c
	)

set	(LIBS_DEPEND
	)

set (LIB Numerical_Flux)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
