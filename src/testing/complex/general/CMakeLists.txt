set	(SOURCE
	 test_complex_const_cast.c
	)

set	(LIBS_DEPEND
	)

set (LIB Test_Complex_General)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
