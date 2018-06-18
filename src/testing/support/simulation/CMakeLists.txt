set	(SOURCE
	 computational_elements/test_support_computational_elements.c
	 computational_elements/test_support_face.c
	 computational_elements/test_support_volume.c

	 solvers/test_support_solve.c
	 solvers/dg/test_support_solve_dg.c
	 solvers/dpg/test_support_solve_dpg.c
	 solvers/opg/test_support_solve_opg.c

	 ../intrusive/test_support_intrusive.c
	)

set	(LIBS_DEPEND
	 Simulation
	 Test_Support_Containers
	)

set (LIB Test_Support_Simulation)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
