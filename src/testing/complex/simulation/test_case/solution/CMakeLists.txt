set	(SOURCE
	 test_complex_solution.c
	 advection/test_complex_solution_advection.c
	 advection/test_complex_solution_advection_default.c
	 advection/peterson/test_complex_solution_peterson.c
	 diffusion/test_complex_solution_diffusion.c
	 diffusion/default_steady/test_complex_solution_diffusion_default_steady.c
	 euler/test_complex_solution_euler.c
	 euler/periodic_vortex/test_complex_solution_periodic_vortex.c
	 euler/supersonic_vortex/test_complex_solution_supersonic_vortex.c
	 euler/free_stream/test_complex_solution_free_stream.c
	 navier_stokes/test_complex_solution_navier_stokes.c
	 navier_stokes/taylor_couette/test_complex_solution_taylor_couette.c
	)

set	(LIBS_DEPEND
	 Test_Complex_Flux
	 Test_Complex_Geometry
	 Test_Complex_Numerical_Flux
	)

set(LIB Test_Complex_Solution)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
