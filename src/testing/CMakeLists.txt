# Library generation
set	(SOURCE
	 test_base.c
	)

set	(LIBS_DEPEND
	 General
	)

set (LIB Test_Base)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})


 include_directories(
 	integration
 	support
 	support/containers
 	support/element
 	support/general
 	support/intrusive
 	support/simulation/computational_elements
 	support/simulation/solvers/dg
 	support/simulation/solvers/dpg
 	unit
 	)

# add_subdirectory(complex)
add_subdirectory(integration)
add_subdirectory(support)
add_subdirectory(unit)
