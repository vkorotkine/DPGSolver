set	(SOURCE
	 multiarray.c
	 matrix.c
	 vector.c
	)

set	(LIBS_DEPEND
	 Allocators
	 m # math
	 ${MKL_LIBRARIES}
	)

add_library(Containers ${SOURCE})
target_link_libraries(Containers PRIVATE ${LIBS_DEPEND})
