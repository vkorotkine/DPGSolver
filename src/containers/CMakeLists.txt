set	(SOURCE
	 multiarray.c
	 multiarray_constructors.c
	 multiarray_math.c
	 multiarray_print.c
	 matrix.c
	 matrix_constructors.c
	 matrix_math.c
	 matrix_print.c
	 vector.c
	 vector_constructors.c
	 vector_math.c
	 vector_print.c
	)

set	(LIBS_DEPEND
	 m # math
	 ${MKL_LIBRARIES}
	 ${GSL_LIBRARIES}
	 General
	 Operators
	)

add_library(Containers STATIC ${SOURCE})
target_link_libraries(Containers PRIVATE ${LIBS_DEPEND})
