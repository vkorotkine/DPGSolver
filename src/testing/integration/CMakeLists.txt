set	(SOURCE
	 test_integration.c
	 test_integration_mesh.c
	 test_integration_fe_init.c
	 test_integration_geometry.c
	 test_integration_euler.c
	 test_integration_linearization.c
	)

set	(LIBS_DEPEND
	 ${MPI_C_LIBRARIES}
	 Geometry
	 Mesh
	 Post_Processing
	 Simulation
	 Test_Base
	 Test_Support
	)

set (LIB Test_Integration)
add_library(${LIB} ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
#print_include_dirs()
