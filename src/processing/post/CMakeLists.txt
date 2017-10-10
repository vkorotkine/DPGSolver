set	(SOURCE
	 element_plotting.c
	 visualization.c
	)

set	(LIBS_DEPEND
	 Computational_Elements
	 Element
	 Operators
	)

add_library(Post_Processing ${SOURCE})
target_link_libraries(Post_Processing PRIVATE ${LIBS_DEPEND})
