set	(SOURCE
	 test_complex_solution.c
	 advection/test_complex_solution_advection.c
	 advection/test_complex_solution_advection_default.c
	 advection/peterson/test_complex_solution_peterson.c
	 euler/test_complex_solution_euler.c
	 euler/periodic_vortex/test_complex_solution_periodic_vortex.c
	 euler/supersonic_vortex/test_complex_solution_supersonic_vortex.c
	)

set	(LIBS_DEPEND
	)

set(LIB Test_Complex_Solution)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
