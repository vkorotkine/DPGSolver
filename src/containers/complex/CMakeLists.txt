set	(SOURCE
	 complex_multiarray.c
	 complex_multiarray_constructors.c
	)

set	(LIBS_DEPEND
	)

set (LIB Complex_Containers)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
