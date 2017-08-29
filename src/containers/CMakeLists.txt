set	(SOURCE
	 multiarray.c
	 matrix.c
	 vector.c
	)

set	(LIBS
	 Allocators
	 m # math
	)

add_library(Containers ${SOURCE})
target_link_libraries(Containers PRIVATE ${LIBS})
