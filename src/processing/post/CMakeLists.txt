set	(SOURCE
	 element_plotting.c
	 visualization.c
	)

set	(LIBS_DEPEND
	 Computational_Elements
	 Element
	 Solvers
	 Operators
	)

add_library(Post_Processing STATIC ${SOURCE})
target_link_libraries(Post_Processing ${LIBS_DEPEND})
