set	(SOURCE
	 test_support_compute_volume_rhs_dg.c
	 test_support_solve_dg.c
	)

set	(LIBS_DEPEND
	 Test_Support_Flux
	 Test_Support_General
	 Test_Support_Operators
	)

set (LIB Test_Support_Solver_DG)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
