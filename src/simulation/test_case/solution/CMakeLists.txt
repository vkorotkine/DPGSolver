set	(SOURCE
	 element_solution.c
	 solution.c
	 advection/solution_advection.c
	 advection/solution_advection_default.c
	 advection/peterson/solution_peterson.c
	 euler/solution_euler.c
	 euler/periodic_vortex/solution_periodic_vortex.c
	 euler/supersonic_vortex/solution_supersonic_vortex.c
	 euler/free_stream/solution_free_stream.c
	)

set	(LIBS_DEPEND
	 Boundary
	 Computational_Elements
	 Compute_Error
	 Containers
	 Element
	 Flux
	 Numerical_Flux
	)

set(LIB Solution)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
