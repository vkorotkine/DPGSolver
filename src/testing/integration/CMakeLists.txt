set	(SOURCE
	 test_integration.c
	)

set	(LIBS_DEPEND
	 Adaptation
	 Simulation
	 Test_Complex_Computational_Elements
	)

set (LIB Test_Integration)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})


set (LIBS_BASE Test_Base Test_Support)

set (EXEC test_integration_mesh)
set (LIBS_DEPEND ${LIBS_BASE} General Mesh Simulation Test_Support_Containers)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "blended_2d_mixed.msh" "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_ParametricMixed2D__ml0__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "straight_2d_quad_periodic.msh" "euler/periodic_vortex/TEST_Euler_PeriodicVortex_QUAD__ml0__p2")

set (EXEC test_integration_fe_init)
set (LIBS_DEPEND ${LIBS_BASE} ${PETSC_LIBRARIES} Simulation Test_Support_Intrusive)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_blended_2d_mixed")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_straight_2d_quad_periodic")

set (EXEC test_integration_geometry)
set (LIBS_DEPEND ${LIBS_BASE} ${PETSC_LIBRARIES} Simulation Test_Support_Containers Test_Support_Intrusive)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_straight_2d_quad_periodic")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_blended_2d_mixed")

set (EXEC test_integration_fluxes)
set (LIBS_DEPEND Test_Base Test_Integration Test_Complex_Test_Case Test_Complex_Solvers)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "flux_advection_1d" "integration/TEST_Advection_Default_1d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "flux_advection_2d" "integration/TEST_Advection_Default_2d__ml0__p2")
#add_test_DPG_w_path(${BIN_PATH_3D} ${EXEC} "flux_advection_3d" "integration/TEST_Advection_Default_3d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "flux_diffusion_1d" "integration/TEST_Diffusion_Default_1d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "flux_diffusion_2d" "integration/TEST_Diffusion_Default_2d__ml0__p2")
#add_test_DPG_w_path(${BIN_PATH_3D} ${EXEC} "flux_diffusion_3d" "integration/TEST_Diffusion_Default_3d__ml0__p2")
#add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "flux_euler_1d" "integration/TEST_Euler_Default_1d__ml0__p0")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "flux_euler_2d" "integration/TEST_Euler_Default_2d__ml0__p0")
#add_test_DPG_w_path(${BIN_PATH_3D} ${EXEC} "flux_euler_3d" "integration/TEST_Euler_Default_3d__ml0__p0")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "flux_navier_stokes_2d" "integration/TEST_NavierStokes_Default_2d__ml0__p0") # TODO: Change this to use DPG such that Hessians are also tested.

set (EXEC test_integration_boundary_conditions)
set (LIBS_DEPEND Test_Base Test_Integration Test_Complex_Test_Case Test_Complex_Solvers)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "boundary_conditions_advection_1d" "integration/TEST_Advection_Default_1d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "boundary_conditions_advection_2d" "integration/TEST_Advection_Default_2d__ml0__p2")
#add_test_DPG_w_path(${BIN_PATH_3D} ${EXEC} "boundary_conditions_advection_3d" "integration/TEST_Advection_Default_3d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "boundary_conditions_diffusion_1d" "integration/TEST_Diffusion_Default_1d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "boundary_conditions_diffusion_2d" "integration/TEST_Diffusion_Default_2d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "boundary_conditions_euler_2d_riemann_slipwall" "integration/boundary_conditions/TEST_Euler_RiSl_2d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "boundary_conditions_euler_2d_supersonic_in/out" "integration/boundary_conditions/TEST_Euler_SiSo_2d__ml0__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "boundary_conditions_navier_stokes_2d_noslip_general/adiabatic" "integration/boundary_conditions/TEST_NavierStokes_NSGA_2d__ml0__p2")

# Potentially add test for equivalence between explicit/implicit computation of grad_coef.

set (EXEC test_integration_conservation)
set (LIBS_DEPEND ${LIBS_BASE} Core Simulation Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/dg/TEST_Advection_Default_DG_Mixed2D__ml2__p2" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/dpg/TEST_Advection_Default_Conservative_DPG_Mixed2D__ml2__p2" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_ParametricMixed2D__ml0__p3" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_Conservative_DPG_BlendedMixed2D__ml1__p2" "petsc_options_gmres_default")

set (EXEC test_integration_linearization)
set (LIBS_DEPEND ${LIBS_BASE} 
     Simulation 
     Test_Integration 
     Test_Complex_Geometry
     Test_Support_General
     Test_Complex_Solution
     Test_Complex_Solvers
     Test_Complex_Test_Case
    )
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "advection/default/dpg/TEST_Advection_Demkowicz_DPGII_LINE__ml2__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/peterson/dg/TEST_Advection_Peterson_TRI__ml2__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/peterson/dpg/TEST_Advection_Peterson_DPG_TRI__ml1__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/dpg/TEST_Advection_Default_Conservative_DPG_Mixed2D__ml2__p2")
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "diffusion/steady/default/dg/TEST_Diffusion_Steady_Default_DG_LINE__ml2__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "diffusion/steady/default/dg/TEST_Diffusion_Steady_Default_DG_TRI__ml2__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "diffusion/steady/default/dg/TEST_Diffusion_Steady_Default_DG_QUAD__ml2__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "diffusion/steady/default/dg/TEST_Diffusion_Steady_Default_DG_Mixed2D__ml2__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_ParametricMixed2D__ml0__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DPG_ParametricMixed2D__ml0__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_BlendedMixed2D__ml0__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_Conservative_DPG_BlendedMixed2D__ml1__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "navier_stokes/steady/taylor_couette/dg/TEST_NavierStokes_TaylorCouette_DG_ParametricMixed2D__ml0__p2")

set (EXEC test_integration_convergence)
set (LIBS_DEPEND ${LIBS_BASE} Core Simulation Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/peterson/dg/TEST_Advection_Peterson_TRI__ml0__p0" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/peterson/dpg/TEST_Advection_Peterson_DPG_TRI__ml0__p0" "petsc_options_cg_ilu1") # Incomplete Cholesky was slower than LU.
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "advection/default/dg/TEST_Advection_Default_DG_LINE" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "advection/default/dpg/TEST_Advection_Demkowicz_DPGII_LINE" "petsc_options_cg_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/dpg/TEST_Advection_Default_DPG_TRI" "petsc_options_cg_ilu1") 
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/dpg/TEST_Advection_Default_DPG_QUAD" "petsc_options_cg_icc")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/dpg/TEST_Advection_Default_DPG_Mixed2D" "petsc_options_cg_icc")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/dpg/TEST_Advection_Default_Conservative_DPG_Mixed2D__ml0" "petsc_options_cg_default")
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "diffusion/steady/default/dg/TEST_Diffusion_Steady_Default_DG_LINE" "petsc_options_cg_ilu1")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "diffusion/steady/default/dg/TEST_Diffusion_Steady_Default_DG_Mixed2D__ml0" "petsc_options_cg_ilu1")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/periodic_vortex/TEST_Euler_PeriodicVortex_QUAD__ml0__p2" "petsc_options_empty")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_ParametricMixed2D" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_Conservative_DPG_BlendedMixed2D__ml1" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DPG_ParametricMixed2D__ml0" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_BlendedMixed2D__ml1" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/joukowski/TEST_Euler_Joukowski_DG_ParametricMixed2D" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DPG_BlendedMixed2D__ml1" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/free_stream/TEST_Euler_FreeStream_DG_StraightTRI2D__ml0" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/free_stream/TEST_Euler_FreeStream_DG_ParametricMixed2D__ml0" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "navier_stokes/steady/taylor_couette/dg/TEST_NavierStokes_TaylorCouette_DG_ParametricTRI" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "navier_stokes/steady/taylor_couette/dg/TEST_NavierStokes_TaylorCouette_DG_ParametricQUAD" "petsc_options_gmres_default")

# Add tests:
# - equivalence operators.
