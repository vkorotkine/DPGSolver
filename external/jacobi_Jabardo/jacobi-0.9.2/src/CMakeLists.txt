set	(SOURCE
	 gauss_quad.c
	 interface.c 
	 jacobi.c
	)

add_library(Jacobi_GSL STATIC ${SOURCE})
target_link_libraries(Jacobi_GSL ${GSL_LIBRARIES})
target_include_directories(Jacobi_GSL PUBLIC .)
