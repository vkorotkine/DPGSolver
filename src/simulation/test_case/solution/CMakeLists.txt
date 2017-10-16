set	(SOURCE
	 element_solution.c
	 solution.c
	 solution_euler.c
	 periodic_vortex/solution_periodic_vortex.c
	)

set	(LIBS_DEPEND
	 Containers
	)

add_library(Solution ${SOURCE})
target_link_libraries(Solution PRIVATE ${LIBS_DEPEND})
