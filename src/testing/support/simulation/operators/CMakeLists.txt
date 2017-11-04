set	(SOURCE
	 test_support_operators.c
	)

set	(LIBS_DEPEND
	 Test_Support_Containers
	)

set (LIB Test_Support_Operators)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
