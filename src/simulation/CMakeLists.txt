add_subdirectory(computational_elements)
add_subdirectory(element)
add_subdirectory(geometry)
add_subdirectory(mesh)
add_subdirectory(solution)

set	(SOURCE
	 simulation.c
	)

set	(LIBS_DEPEND
	 Computational_Elements
	)

add_library(Simulation ${SOURCE})
target_link_libraries(Simulation PRIVATE ${LIBS_DEPEND})
