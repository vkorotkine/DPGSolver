set	(SOURCE
	 element.c
	 element_operators.c
	 bases.c
	 cubature.c
	)

set	(LIBS_DEPEND
	 ${GSL_LIBRARIES}
	 Containers
	 Intrusive
	 General
	 Jacobi_GSL
	)

add_library(Element STATIC ${SOURCE})
target_link_libraries(Element PRIVATE ${LIBS_DEPEND})
