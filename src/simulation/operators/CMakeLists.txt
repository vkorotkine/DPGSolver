set	(SOURCE
	 multiarray_operator.cpp
	 operator.cpp
	)

set	(LIBS_DEPEND
	 Containers
	)

add_library(Operators STATIC ${SOURCE})
target_link_libraries(Operators PRIVATE ${LIBS_DEPEND})
