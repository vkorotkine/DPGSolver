set	(SOURCE
	 test_support_computational_elements.c
	 test_support_face.c
	 test_support_volume.c
	)

set	(LIBS_DEPEND
	 Test_Support_Containers
	)

set (LIB Test_Support_Computational_Elements)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
