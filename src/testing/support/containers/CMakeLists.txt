set	(SOURCE
	 test_support_matrix.c
	 test_support_multiarray.c
	 test_support_vector.c
	)

set	(LIBS_DEPEND
	 Containers
	 Test_Support
	)

set (LIB Test_Support_Containers)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
