set	(SOURCE
	 test_support_bases.cpp
	 test_support_nodes.cpp
	)

set	(LIBS_DEPEND
	 Element
	)

set (LIB Test_Support_Element)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
