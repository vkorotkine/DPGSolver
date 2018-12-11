add_subdirectory(mesh)
add_subdirectory(operators)

set	(SOURCE
	 simulation.cpp

	 computational_elements/computational_elements.cpp
	 computational_elements/face.cpp
	 computational_elements/volume.cpp

	 geometry/element_geometry.cpp
	 geometry/geometry.cpp
	 geometry/geometry_blended.cpp
	 geometry/geometry_parametric.cpp
	 geometry/geometry_surface.cpp
	 geometry/geometry_normals.cpp

	 solvers/compute_rlhs.cpp
	 solvers/compute_face_rlhs.cpp
	 solvers/compute_volume_rlhs.cpp
	 solvers/element_solver.cpp
	 solvers/face_solver.cpp
	 solvers/solve.cpp
	 solvers/solve_explicit.cpp
	 solvers/solve_implicit.cpp
	 solvers/volume_solver.cpp

	 solvers/adaptation/adaptation.cpp
	 solvers/adaptation/face_solver_adaptive.cpp
	 solvers/adaptation/element_adaptation.cpp
	 solvers/adaptation/volume_solver_adaptive.cpp

	 solvers/dg/compute_face_rlhs_dg.cpp
	 solvers/dg/compute_grad_coef_dg.cpp
	 solvers/dg/compute_source_rlhs_dg.cpp
	 solvers/dg/compute_volume_rlhs_dg.cpp
	 solvers/dg/element_solver_dg.cpp
	 solvers/dg/face_solver_dg.cpp
	 solvers/dg/solve_dg.cpp
	 solvers/dg/volume_solver_dg.cpp

	 solvers/dpg/compute_all_rlhs_dpg.cpp
	 solvers/dpg/face_solver_dpg.cpp
	 solvers/dpg/element_solver_dpg.cpp
	 solvers/dpg/solve_dpg.cpp
	 solvers/dpg/volume_solver_dpg.cpp

	 solvers/opg/compute_volume_rlhs_opg.cpp
	 solvers/opg/compute_face_rlhs_opg.cpp
	 solvers/opg/face_solver_opg.cpp
	 solvers/opg/element_solver_opg.cpp
	 solvers/opg/penalty_opg.cpp
	 solvers/opg/penalty_opg_advection.cpp
	 solvers/opg/solve_opg.cpp
	 solvers/opg/volume_solver_opg.cpp

	 solvers/numerical_flux/numerical_flux.cpp
	 solvers/numerical_flux/numerical_flux_advection.cpp
	 solvers/numerical_flux/numerical_flux_diffusion.cpp
	 solvers/numerical_flux/numerical_flux_euler.cpp
	 solvers/numerical_flux/numerical_flux_navier_stokes.cpp
	 solvers/numerical_flux/numerical_flux_burgers_inviscid.cpp

	 test_case/test_case.cpp
	 test_case/test_case_c.cpp

	 test_case/boundary/boundary.cpp
	 test_case/boundary/boundary_advection.cpp
	 test_case/boundary/boundary_diffusion.cpp
	 test_case/boundary/boundary_euler.cpp
	 test_case/boundary/boundary_navier_stokes.cpp

	 test_case/compute_error/compute_error.cpp
	 test_case/compute_error/compute_error_advection.cpp
	 test_case/compute_error/compute_error_diffusion.cpp
	 test_case/compute_error/compute_error_euler.cpp
	 test_case/compute_error/compute_error_navier_stokes.cpp
	 test_case/compute_error/compute_error_burgers_inviscid.cpp

	 test_case/flux/flux.cpp
	 test_case/flux/flux_advection.cpp
	 test_case/flux/flux_diffusion.cpp
	 test_case/flux/flux_euler.cpp
	 test_case/flux/flux_navier_stokes.cpp
	 test_case/flux/flux_burgers_inviscid.cpp

	 test_case/solution/element_solution.cpp
	 test_case/solution/solution.cpp
	 test_case/solution/advection/solution_advection.cpp
	 test_case/solution/advection/solution_advection_default.cpp
	 test_case/solution/advection/free_stream_advection/solution_free_stream_advection.cpp
	 test_case/solution/advection/peterson/solution_peterson.cpp
	 test_case/solution/advection/vortex_advection/solution_vortex_advection.cpp
	 test_case/solution/advection/hyperbolic_tan/solution_hyperbolic_tan.cpp
	 test_case/solution/diffusion/solution_diffusion.cpp
	 test_case/solution/diffusion/default_steady/solution_diffusion_default_steady.cpp
	 test_case/solution/euler/solution_euler.cpp
	 test_case/solution/euler/periodic_vortex/solution_periodic_vortex.cpp
	 test_case/solution/euler/supersonic_vortex/solution_supersonic_vortex.cpp
	 test_case/solution/euler/free_stream/solution_free_stream.cpp
	 test_case/solution/navier_stokes/solution_navier_stokes.cpp
	 test_case/solution/navier_stokes/taylor_couette/solution_taylor_couette.cpp
	 test_case/solution/burgers_inviscid/solution_burgers_inviscid.cpp
	 test_case/solution/burgers_inviscid/trigonometric/solution_trigonometric.cpp
	 test_case/solution/restart/approximate_nearest_neighbor.cpp
	 test_case/solution/restart/inverse_mapping.cpp
	 test_case/solution/restart/restart.cpp
	 test_case/solution/restart/restart_writers.cpp

	 visualization/element_plotting.cpp
	 visualization/visualization.cpp
	)

set	(LIBS_DEPEND
	 ${MPI_C_LIBRARIES}
	 ${PETSC_LIBRARIES}
	 Containers
	 Intrusive
	 Mesh
	)

set(LIB Simulation)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
