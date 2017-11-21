add_subdirectory(computational_elements)
add_subdirectory(geometry)
add_subdirectory(operators)
add_subdirectory(mesh)
add_subdirectory(solvers)
add_subdirectory(test_case)

set	(SOURCE
	 simulation.c
	)

set	(LIBS_DEPEND
	 ${MPI_C_LIBRARIES}
	 Computational_Elements
	 Mesh
	 Test_Case
	)

set(LIB Simulation)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
