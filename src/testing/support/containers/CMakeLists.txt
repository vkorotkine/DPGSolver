set	(SOURCE
	 test_support_matrix.cpp
	 test_support_multiarray.cpp
	 test_support_vector.cpp
	)

set	(LIBS_DEPEND
	 Containers
	 Test_Support
	)

set (LIB Test_Support_Containers)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
