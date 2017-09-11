set	(SOURCE
	 geometry.c
	 geometry_element.c
	)

set	(LIBS_DEPEND
	 Computational_Elements
	 Containers
	 Intrusive
	)

add_library(Geometry STATIC ${SOURCE})
target_link_libraries(Geometry PRIVATE ${LIBS_DEPEND})
