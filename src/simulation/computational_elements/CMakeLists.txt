set	(SOURCE
	 computational_elements.c
	 face.c
	 face_solver.c
	 volume.c
	 volume_solver.c
	)

set	(LIBS_DEPEND
	 Containers
	 Element
	 General
	 Geometry
	 Intrusive
	 Post_Processing
	 Solution
	)

set(LIB Computational_Elements)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
