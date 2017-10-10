set	(SOURCE
	 test_unit.c
	 test_unit_bases.c
	 test_unit_containers.c
	 test_unit_nodes.c
	 test_unit_operators_tp.c
	)

set	(LIBS_DEPEND
	 Test_Base
	 Test_Support
	)

set (LIB Test_Unit)
add_library(${LIB} ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
