set	(SOURCE
	 test_complex_const_cast.c
	 test_complex_math_functions.c
	)

set	(LIBS_DEPEND
	)

set (LIB Test_Complex_General)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
