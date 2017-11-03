add_subdirectory(containers)
add_subdirectory(element)
add_subdirectory(general)
add_subdirectory(intrusive)
add_subdirectory(simulation)

set	(SOURCE
	 test_support.c
	)

### TODO: Needs clean-up for the library dependencies.
set	(LIBS_DEPEND
	 ${GSL_LIBRARIES}
	 Complex_Containers
	 Computational_Elements
	 Containers
	 Element
	 General
	 Intrusive
	 Test_Support_Containers
	 Test_Support_Element
	 Test_Support_Flux
	 Test_Support_General
	 Test_Support_Intrusive
	 Test_Support_Computational_Elements
	 Test_Support_Solver_DG
	)

add_library(Test_Support ${SOURCE})
target_link_libraries(Test_Support PRIVATE ${LIBS_DEPEND})
