set	(SOURCE
	 computational_elements/test_support_computational_elements.cpp
	 computational_elements/test_support_face.cpp
	 computational_elements/test_support_volume.cpp

	 solvers/test_support_solve.cpp
	 solvers/dg/test_support_solve_dg.cpp
	 solvers/dpg/test_support_solve_dpg.cpp
	 solvers/opg/test_support_solve_opg.cpp

	 ../intrusive/test_support_intrusive.cpp
	)

set	(LIBS_DEPEND
	 Simulation
	 Test_Support_Containers
	)

set (LIB Test_Support_Simulation)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
