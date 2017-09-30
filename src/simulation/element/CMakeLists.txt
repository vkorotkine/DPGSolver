set	(SOURCE
	 element.c
	 element_operators.c
	 element_operators_tp.c
	 bases.c
	 cubature.c
	 cubature_operators.c
	)

message("${GSL_LIBRARIES}")
set	(LIBS_DEPEND
	 ${GSL_LIBRARIES}
	 ${MKL_LIBRARIES} # Likely to be removed
	 Containers
	 Intrusive
	 General
	 Jacobi_GSL
	)

add_library(Element STATIC ${SOURCE})
target_link_libraries(Element PRIVATE ${LIBS_DEPEND})
