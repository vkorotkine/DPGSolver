set	(SOURCE
	 test_complex_boundary.c
	 test_complex_boundary_advection.c
	)

set	(LIBS_DEPEND
	)

set (LIB Test_Complex_Boundary)
add_library(${LIB} ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
