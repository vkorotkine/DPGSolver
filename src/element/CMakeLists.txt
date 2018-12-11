set	(SOURCE
	 bases.cpp
	 element.cpp
	 element_operators.cpp
	 element_operators_tp.cpp
	 nodes.cpp
	 nodes_correspondence.cpp
	 nodes_operators.cpp
	 nodes_plotting.cpp
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
