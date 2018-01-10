set	(SOURCE
	 test_support_intrusive.c
	)

set	(LIBS_DEPEND
	 Computational_Elements
	 Test_Support_Computational_Elements
	)

set (LIB Test_Support_Intrusive)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
