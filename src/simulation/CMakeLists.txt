include_directories(
	mesh
	)

add_subdirectory(mesh)

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
	 face.c
	)

set (LIBS_DEPEND
	 Allocators
	 Containers
	 Element
	 General
	 Intrusive
	)

add_library(Simulation ${SOURCE})
target_link_libraries(Simulation PRIVATE ${LIBS_DEPEND})
