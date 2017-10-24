set	(SOURCE
	 compute_error.c
	)

set	(LIBS_DEPEND
	)

add_library(Compute_Error ${SOURCE})
target_link_libraries(Compute_Error PRIVATE ${LIBS_DEPEND})
