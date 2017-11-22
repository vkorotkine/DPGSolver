set	(SOURCE
	 test_integration.c
	)

set	(LIBS_DEPEND
	 Adaptation
	 Simulation
	)

set (LIB Test_Integration)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})


set (LIBS_BASE Test_Base Test_Support)

set (EXEC test_integration_mesh)
set (LIBS_DEPEND ${LIBS_BASE} General Mesh Test_Support_Containers)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} PRIVATE ${LIBS_DEPEND})
add_test_DPG(${EXEC} "curved_2d_mixed.msh")
add_test_DPG(${EXEC} "straight_2d_quad_periodic.msh")

set (EXEC test_integration_fe_init)
set (LIBS_DEPEND ${LIBS_BASE} ${PETSC_LIBRARIES} Simulation Test_Support_Intrusive)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} PRIVATE ${LIBS_DEPEND})
add_test_DPG(${EXEC} "extern_mesh/TEST_curved_2d_mixed")
add_test_DPG(${EXEC} "extern_mesh/TEST_straight_2d_quad_periodic")

set (EXEC test_integration_geometry)
set (LIBS_DEPEND ${LIBS_BASE} ${PETSC_LIBRARIES} Simulation Test_Support_Containers Test_Support_Intrusive)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} PRIVATE ${LIBS_DEPEND})
add_test_DPG(${EXEC} "extern_mesh/TEST_straight_2d_quad_periodic")
add_test_DPG(${EXEC} "extern_mesh/TEST_curved_2d_mixed")

set (EXEC test_integration_linearization)
set (LIBS_DEPEND ${LIBS_BASE} Simulation Test_Complex_Solvers Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} PRIVATE ${LIBS_DEPEND})
add_test_DPG(${EXEC} "advection/TEST_Advection_Demkowicz_DPGII_LINE__ml0__p1")
add_test_DPG(${EXEC} "advection/TEST_Advection_Peterson_TRI__ml0__p1")
add_test_DPG(${EXEC} "advection/TEST_Advection_Peterson_DPG_TRI__ml0")

set (EXEC test_integration_convergence)
set (LIBS_DEPEND ${LIBS_BASE} Simulation Test_Integration)
add_executable(${EXEC} ${EXEC}.c)
target_link_libraries(${EXEC} PRIVATE ${LIBS_DEPEND})
add_test_DPG(${EXEC} "advection/TEST_Advection_Peterson_TRI__ml0__p1")
add_test_DPG(${EXEC} "advection/TEST_Advection_Demkowicz_DPGII_LINE__ml0__p1")
add_test_DPG(${EXEC} "euler/TEST_Euler_PeriodicVortex_QUAD__ml0__p2")

# Add tests:
# - equivalence real/complex.
# - equivalence operators.
