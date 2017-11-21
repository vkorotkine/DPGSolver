add_subdirectory(dg)
add_subdirectory(dpg)
add_subdirectory(numerical_flux)

set	(SOURCE
	 test_complex_compute_face_rhs.c
	 test_complex_compute_volume_rhs.c
	)

set	(LIBS_DEPEND
#	 Test_Complex_Numerical_Flux
	 Test_Complex_Solver_DG
	 Test_Complex_Solver_DPG
	)

set (LIB Test_Complex_Solvers)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
