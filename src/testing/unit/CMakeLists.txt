set	(SOURCE
	 test_unit.c
	)

set	(LIBS_DEPEND
	 Mesh
	)

set (LIB Test_Unit)
add_library(${LIB} ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
