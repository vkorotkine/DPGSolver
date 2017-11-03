set	(SOURCE
	 test_support_math_functions.c
	)

set	(LIBS_DEPEND
	)

set (LIB Test_Support_General)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
