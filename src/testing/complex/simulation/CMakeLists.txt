add_subdirectory(operators)

set	(SOURCE
	 computational_elements/test_complex_computational_elements.c

	 geometry/test_complex_geometry.c
	 geometry/test_complex_geometry_blended.c
	 geometry/test_complex_geometry_parametric.c
	 geometry/test_complex_geometry_surface.c

	 solvers/test_complex_compute_face_rhs.c
	 solvers/test_complex_compute_volume_rhs.c
	 solvers/test_complex_face_solver.c
	 solvers/test_complex_solve.c
	 solvers/test_complex_volume_solver.c

	 solvers/dg/test_complex_compute_grad_coef_dg.c
	 solvers/dg/test_complex_compute_face_rhs_dg.c
	 solvers/dg/test_complex_compute_volume_rhs_dg.c
	 solvers/dg/test_complex_solve_dg.c
	 solvers/dg/test_complex_face_solver_dg.c
	 solvers/dg/test_complex_volume_solver_dg.c

	 solvers/dpg/test_complex_compute_all_rhs_dpg.c
	 solvers/dpg/test_complex_solve_dpg.c
	 solvers/dpg/test_complex_face_solver_dpg.c
	 solvers/dpg/test_complex_volume_solver_dpg.c

	 solvers/numerical_flux/test_complex_numerical_flux.c
	 solvers/numerical_flux/test_complex_numerical_flux_advection.c
	 solvers/numerical_flux/test_complex_numerical_flux_diffusion.c
	 solvers/numerical_flux/test_complex_numerical_flux_euler.c
	 solvers/numerical_flux/test_complex_numerical_flux_navier_stokes.c

	 test_case/test_complex_test_case.c

	 test_case/boundary/test_complex_boundary.c
	 test_case/boundary/test_complex_boundary_advection.c
	 test_case/boundary/test_complex_boundary_diffusion.c
	 test_case/boundary/test_complex_boundary_euler.c
	 test_case/boundary/test_complex_boundary_navier_stokes.c

	 test_case/flux/test_complex_flux.c
	 test_case/flux/test_complex_flux_advection.c
	 test_case/flux/test_complex_flux_diffusion.c
	 test_case/flux/test_complex_flux_euler.c
	 test_case/flux/test_complex_flux_navier_stokes.c

	 test_case/solution/test_complex_solution.c
	 test_case/solution/advection/test_complex_solution_advection.c
	 test_case/solution/advection/test_complex_solution_advection_default.c
	 test_case/solution/advection/peterson/test_complex_solution_peterson.c
	 test_case/solution/advection/free_stream_advection/test_complex_solution_free_stream_advection.c
	 test_case/solution/advection/vortex_advection/test_complex_solution_vortex_advection.c
	 test_case/solution/diffusion/test_complex_solution_diffusion.c
	 test_case/solution/diffusion/default_steady/test_complex_solution_diffusion_default_steady.c
	 test_case/solution/euler/test_complex_solution_euler.c
	 test_case/solution/euler/periodic_vortex/test_complex_solution_periodic_vortex.c
	 test_case/solution/euler/supersonic_vortex/test_complex_solution_supersonic_vortex.c
	 test_case/solution/euler/free_stream/test_complex_solution_free_stream.c
	 test_case/solution/navier_stokes/test_complex_solution_navier_stokes.c
	 test_case/solution/navier_stokes/taylor_couette/test_complex_solution_taylor_couette.c
	 test_case/solution/restart/test_complex_restart.c
	)

set	(LIBS_DEPEND
	 ${MPI_C_LIBRARIES}
	 ${PETSC_LIBRARIES}
	 Containers
	 Intrusive
	 Mesh
	 Test_Complex_Containers
	 Test_Support_Containers
	 Test_Support_Simulation
	)

set(LIB Test_Complex_Simulation)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
