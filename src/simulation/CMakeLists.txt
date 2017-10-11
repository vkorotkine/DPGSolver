add_subdirectory(computational_elements)
add_subdirectory(geometry)
add_subdirectory(operators)
add_subdirectory(mesh)
add_subdirectory(solution)

set	(SOURCE
	 simulation.c
	)

set	(LIBS_DEPEND
	 ${MPI_C_LIBRARIES}
	 Computational_Elements
	 Mesh
	)

add_library(Simulation ${SOURCE})
target_link_libraries(Simulation PRIVATE ${LIBS_DEPEND})
