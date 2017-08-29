set	(SOURCE
	test_integration_euler.c
	test_integration_linearization.c
	)

set	(LIBS
	 Simulation
	 Mesh
	)

add_library(Test_Integration ${SOURCE})
target_link_libraries(Test_Integration PRIVATE ${LIBS})
#print_include_dirs()
