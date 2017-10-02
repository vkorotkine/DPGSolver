set	(SOURCE
	 multiarray_operator.c
	 operator.c
	)

set	(LIBS_DEPEND
	)

add_library(Operators ${SOURCE})
target_link_libraries(Operators PRIVATE ${LIBS_DEPEND})
