set	(SOURCE
	 computational_elements.c
	 volume.c
	 solver_volume.c
	 solver_face.c
	 face.c
	)

set	(LIBS_DEPEND
	 Containers
	 Element
	 General
	 Geometry
	 Intrusive
	)

set(LIB Computational_Elements)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
