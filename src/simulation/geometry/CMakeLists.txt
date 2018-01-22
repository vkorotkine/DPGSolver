set	(SOURCE
	 element_geometry.c
	 geometry.c
	 geometry_blended.c
	 geometry_parametric.c
	)

set	(LIBS_DEPEND
	 Computational_Elements
	 Containers
	 Intrusive
	)

add_library(Geometry STATIC ${SOURCE})
target_link_libraries(Geometry ${LIBS_DEPEND})
