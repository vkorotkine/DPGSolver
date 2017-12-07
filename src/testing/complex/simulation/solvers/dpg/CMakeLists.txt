set	(SOURCE
	 test_complex_compute_all_rhs_dpg.c
	 test_complex_solve_dpg.c
	 test_complex_face_solver_dpg.c
	 test_complex_volume_solver_dpg.c
	)

set	(LIBS_DEPEND
	 Complex_Containers
	 Test_Complex_Boundary
	 Test_Complex_Containers
	 Test_Complex_Flux
	 Test_Complex_Numerical_Flux
	 Test_Complex_Operators
	 Test_Complex_Solvers
	 Test_Complex_Test_Case
	 Test_Support_Containers
	 Test_Support_General
	)

set (LIB Test_Complex_Solver_DPG)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
