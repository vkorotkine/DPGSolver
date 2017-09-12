set	(SOURCE
	 const_cast.c
	 file_processing.c
	 math_functions.c
	)

add_library(General ${SOURCE})
target_link_libraries(General ${GSL_LIBRARIES})
