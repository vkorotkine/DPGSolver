set	(SOURCE
	 test_unit.c
	 test_unit_containers.c
	)

set	(LIBS_DEPEND
	)

set (LIB Test_Unit)
add_library(${LIB} ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
