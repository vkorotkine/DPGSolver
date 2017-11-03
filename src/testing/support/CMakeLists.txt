set	(SOURCE
	 solvers/dg/test_support_solve_dg.c
	 solvers/dg/test_support_compute_volume_rhs_dg.c
	 test_case/flux/test_support_flux.c
	 test_case/flux/test_support_flux_advection.c
	 test_support.c
	 test_support_bases.c
	 test_support_face.c
	 test_support_intrusive.c
	 test_support_math_functions.c
	 test_support_matrix.c
	 test_support_multiarray.c
	 test_support_nodes.c
	 test_support_vector.c
	 test_support_volume.c
	)

set	(LIBS_DEPEND
	 ${GSL_LIBRARIES}
	 Complex_Containers
	 Computational_Elements
	 Containers
	 Element
	 General
	 Intrusive
	)

add_library(Test_Support ${SOURCE})
target_link_libraries(Test_Support ${LIBS_DEPEND})
