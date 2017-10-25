set	(SOURCE
	 compute_error.c
	 compute_error_euler.c
	 element_error.c
	)

set	(LIBS_DEPEND
	 Computational_Elements
	 Containers
	 Element
	)

add_library(Compute_Error STATIC ${SOURCE})
target_link_libraries(Compute_Error PRIVATE ${LIBS_DEPEND})
