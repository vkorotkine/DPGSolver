set	(SOURCE
	 computational_elements.c
	 face.c
	 volume.c
	)

set	(LIBS_DEPEND
	 Compute_Error
	 Containers
	 Element
	 Geometry
	 General
	 Intrusive
	 Post_Processing
	 Solution
	 Solvers
	)

set(LIB Computational_Elements)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
