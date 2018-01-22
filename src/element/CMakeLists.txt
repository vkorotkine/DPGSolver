set	(SOURCE
	 bases.c
	 element.c
	 element_operators.c
	 element_operators_tp.c
	 nodes.c
	 nodes_correspondence.c
	 nodes_operators.c
	 nodes_plotting.c
	)

set	(LIBS_DEPEND
	 ${GSL_LIBRARIES}
	 ${MKL_LIBRARIES} # Likely to be removed
	 Containers
	 Intrusive
	 General
	 Jacobi_GSL
	 Operators
	)

add_library(Element STATIC ${SOURCE})
target_link_libraries(Element ${LIBS_DEPEND})
