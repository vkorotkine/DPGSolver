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
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})


set (LIBS_BASE Test_Base Test_Support)

set (EXEC test_integration_mesh)
set (LIBS_DEPEND ${LIBS_BASE} General Mesh Test_Support_Containers)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} PRIVATE ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "curved_2d_mixed.msh")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "straight_2d_quad_periodic.msh")

set (EXEC test_integration_fe_init)
set (LIBS_DEPEND ${LIBS_BASE} ${PETSC_LIBRARIES} Simulation Test_Support_Intrusive)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} PRIVATE ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_curved_2d_mixed")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_straight_2d_quad_periodic")

set (EXEC test_integration_geometry)
set (LIBS_DEPEND ${LIBS_BASE} ${PETSC_LIBRARIES} Simulation Test_Support_Containers Test_Support_Intrusive)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} PRIVATE ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_straight_2d_quad_periodic")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "extern_mesh/TEST_curved_2d_mixed")

set (EXEC test_integration_fluxes)
set (LIBS_DEPEND Test_Base Test_Integration Test_Complex_Test_Case)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} PRIVATE ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "flux_advection_1d" "integration/fluxes/TEST_Advection_Default_1d__ml0__p0")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "flux_advection_2d" "integration/fluxes/TEST_Advection_Default_2d__ml0__p0")
add_test_DPG_w_path(${BIN_PATH_3D} ${EXEC} "flux_advection_3d" "integration/fluxes/TEST_Advection_Default_3d__ml0__p0")
#add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "flux_euler_1d" "integration/fluxes/TEST_Euler_Default_1d__ml0__p0")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "flux_euler_2d" "integration/fluxes/TEST_Euler_Default_2d__ml0__p0")
#add_test_DPG_w_path(${BIN_PATH_3D} ${EXEC} "flux_euler_3d" "integration/fluxes/TEST_Euler_Default_3d__ml0__p0")

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
target_link_libraries(${EXEC} PRIVATE ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "advection/TEST_Advection_Demkowicz_DPGII_LINE__ml2__p2")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/peterson/dg/TEST_Advection_Peterson_TRI__ml2__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/peterson/dpg/TEST_Advection_Peterson_DPG_TRI__ml1__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_ParametricMixed2D__ml0__p3")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DPG_ParametricMixed2D__ml0__p3")

set (EXEC test_integration_convergence)
set (LIBS_DEPEND ${LIBS_BASE} Core Simulation Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} PRIVATE ${LIBS_DEPEND})
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/peterson/dg/TEST_Advection_Peterson_TRI__ml0__p0" "petsc_options_gmres_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/peterson/dpg/TEST_Advection_Peterson_DPG_TRI__ml0__p0" "petsc_options_cg_ilu1") # Incomplete Cholesky was slower than LU.
add_test_DPG_w_path(${BIN_PATH_1D} ${EXEC} "advection/TEST_Advection_Demkowicz_DPGII_LINE__ml0__p1" "petsc_options_cg_default")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "advection/default/TEST_Advection_Default_DPG_TRI__ml0" "petsc_options_cg_ilu1") # TODO Replace with MIXED2D.
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/periodic_vortex/TEST_Euler_PeriodicVortex_QUAD__ml0__p2" "petsc_options_empty")
add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_ParametricMixed2D__ml0" "petsc_options_gmres_default")
#add_test_DPG_w_path(${BIN_PATH_2D} ${EXEC} "euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DPG_ParametricMixed2D__ml0" "petsc_options_empty")

# Add tests:
# - equivalence operators.
