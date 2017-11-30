set	(SOURCE
	 const_cast.c
	 file_processing.c
	 math_functions.c
	)

set	(LIBS_DEPEND
	 ${GSL_LIBRARIES}
	 Jacobi_GSL
	)

add_library(General STATIC ${SOURCE})
target_link_libraries(General PRIVATE ${LIBS_DEPEND})
