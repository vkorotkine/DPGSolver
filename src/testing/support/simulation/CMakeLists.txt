set	(SOURCE
	 computational_elements/test_support_computational_elements.c
	 computational_elements/test_support_face.c
	 computational_elements/test_support_volume.c

	 ../intrusive/test_support_intrusive.c
	)

set	(LIBS_DEPEND
	 Simulation
	 Test_Support_Containers
	)

set (LIB Test_Support_Simulation)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
