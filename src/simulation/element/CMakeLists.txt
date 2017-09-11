set	(SOURCE
	 element.c
	 cubature.c
	)

set	(LIBS_DEPEND
	 Containers
	 Intrusive
	 General
	)

add_library(Element ${SOURCE})
target_link_libraries(Element PRIVATE ${LIBS_DEPEND})
