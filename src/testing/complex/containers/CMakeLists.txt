set	(SOURCE
	 complex_matrix.c
	 complex_matrix_constructors.c
	 complex_matrix_math.c
	 complex_matrix_print.c
	 complex_multiarray.c
	 complex_multiarray_constructors.c
	 complex_multiarray_math.c
	 complex_multiarray_print.c
	 complex_vector.c
	 complex_vector_constructors.c
	 complex_vector_math.c
	 complex_vector_print.c
	)

set	(LIBS_DEPEND
	 Containers
	)

set (LIB Test_Complex_Containers)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})

set (LIB Complex_Containers)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
