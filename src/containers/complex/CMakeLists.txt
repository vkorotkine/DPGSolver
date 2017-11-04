set	(SOURCE
	 complex_multiarray_minimal_constructors.c
	)

set	(LIBS_DEPEND
	)

set (LIB Complex_Containers)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
