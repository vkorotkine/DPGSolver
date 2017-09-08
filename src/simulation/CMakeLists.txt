add_subdirectory(mesh)
add_subdirectory(geometry)
add_subdirectory(solution)

set (SOURCE
	 element.c
	)

set (LIBS_DEPEND
	 Containers
	 Intrusive
	 General
	)

add_library(Element ${SOURCE})
target_link_libraries(Element PRIVATE ${LIBS_DEPEND})


set	(SOURCE
	 simulation.c
	 volume.c
	 solver_volume.c
	 face.c
	)

set (LIBS_DEPEND
	 Containers
	 Element
	 General
	 Geometry
	 Intrusive
	 Solution
	)

add_library(Simulation ${SOURCE})
target_link_libraries(Simulation PRIVATE ${LIBS_DEPEND})
