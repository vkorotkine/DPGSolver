set	(SOURCE
	 geometry.c
	 geometry_element.c
	)

set	(LIBS_DEPEND
	 Containers
	)

add_library(Geometry ${SOURCE})
target_link_libraries(Geometry PRIVATE ${LIBS_DEPEND})
