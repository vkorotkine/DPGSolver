set	(SOURCE
	 multiarray.cpp
	 multiarray_constructors.cpp
	 multiarray_math.cpp
	 multiarray_print.cpp
	 matrix.cpp
	 matrix_constructors.cpp
	 matrix_math.cpp
	 matrix_print.cpp
	 vector.cpp
	 vector_constructors.cpp
	 vector_math.cpp
	 vector_print.cpp
	)

set	(LIBS_DEPEND
	 m # math
	 ${MKL_LIBRARIES}
	 ${GSL_LIBRARIES}
	 General
	 Operators
	)

add_library(Containers STATIC ${SOURCE})
target_link_libraries(Containers ${LIBS_DEPEND})
