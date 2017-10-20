set	(SOURCE
	 numerical_flux.c
	)

set	(LIBS_DEPEND
	)

set (LIB Numerical_Flux)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
