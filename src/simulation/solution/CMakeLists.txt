set	(SOURCE
	 solution.c
	)

set	(LIBS_DEPEND
	 Containers
	)

add_library(Solution ${SOURCE})
target_link_libraries(Solution PRIVATE ${LIBS_DEPEND})
