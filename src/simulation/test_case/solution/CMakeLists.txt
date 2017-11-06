set	(SOURCE
	 element_solution.c
	 solution.c
	 advection/solution_advection.c
	 advection/solution_advection_default.c
	 advection/peterson/solution_peterson.c
	 euler/solution_euler.c
	 euler/periodic_vortex/solution_periodic_vortex.c
#	 supersonic_vortex/solution_supersonic_vortex.c
	)

set	(LIBS_DEPEND
	 Boundary
	 Computational_Elements
	 Compute_Error
	 Containers
	 Element
	)

add_library(Solution STATIC ${SOURCE})
target_link_libraries(Solution PRIVATE ${LIBS_DEPEND})
