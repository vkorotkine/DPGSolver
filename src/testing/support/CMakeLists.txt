set	(SOURCE
	 test_support.c
	 test_support_vector.c
	 test_support_matrix.c
	 test_support_multiarray.c
	 test_support_intrusive.c
	 test_support_volume.c
	 test_support_face.c
	 test_support_cubature.c
	)

set	(LIBS_DEPEND
	 Containers
	 Element
	 General
	 Intrusive
	)

add_library(Test_Support ${SOURCE})
target_link_libraries(Test_Support ${LIBS_DEPEND})
