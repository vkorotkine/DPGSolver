include_directories(
	mesh
	)

add_subdirectory(mesh)

set(SOURCE simulation.c
           element.c
           volume.c
           face.c
   )

add_library(Simulation ${SOURCE})
