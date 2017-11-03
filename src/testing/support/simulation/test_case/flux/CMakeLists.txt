set	(SOURCE
	 test_support_flux.c
	 test_support_flux_advection.c
	)

set	(LIBS_DEPEND
	)

set (LIB Test_Support_Flux)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
