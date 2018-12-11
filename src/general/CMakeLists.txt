set	(SOURCE
	 const_cast.cpp
	 file_processing.cpp
	 file_processing_conversions.cpp
	 math_functions.cpp
	)

set	(LIBS_DEPEND
	 ${GSL_LIBRARIES}
	 Jacobi_GSL
	)

add_library(General STATIC ${SOURCE})
target_link_libraries(General ${LIBS_DEPEND})
