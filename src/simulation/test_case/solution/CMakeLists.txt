set	(SOURCE
	 element_solution.c
	 solution.c
	 solution_euler.c
	 periodic_vortex/solution_periodic_vortex.c
#	 supersonic_vortex/solution_supersonic_vortex.c
	)

set	(LIBS_DEPEND
	 Computational_Elements
	 Compute_Error
	 Containers
	 Element
	)

add_library(Solution STATIC ${SOURCE})
target_link_libraries(Solution PRIVATE ${LIBS_DEPEND})
