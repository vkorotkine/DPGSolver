set	(SOURCE
	 test_complex_geometry.c
	 test_complex_geometry_blended.c
	 test_complex_geometry_parametric.c
	)

set	(LIBS_DEPEND
	)

set(LIB Test_Complex_Geometry)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
