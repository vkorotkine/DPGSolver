include_directories(
	mesh
	)

add_subdirectory(mesh)

set	(SOURCE
	 simulation.c
	 element.c
	 volume.c
	 face.c
	)

set (LIBS
	 Allocators
	 Containers
	 General
	 Intrusive
	)

add_library(Simulation ${SOURCE})
target_link_libraries(Simulation PRIVATE ${LIBS})
