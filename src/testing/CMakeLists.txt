# Library generation
set	(SOURCE
	 test_base.c
	)

set	(LIBS_DEPEND
	 General
	)

set (LIB Test_Base)
add_library(${LIB} ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})


include_directories(
	support
	support/containers
	support/element
	support/general
	support/intrusive
	support/simulation/computational_elements
	support/simulation/operators
	support/simulation/solvers/dg
	support/simulation/test_case/flux
	unit
	integration
	)

add_subdirectory(support)
add_subdirectory(unit)
add_subdirectory(integration)

# Executable generation
set	(LIBS_DEPEND
	 Test_Base
	 Test_Support
	 Test_Unit 
	 Test_Integration
	)

add_executable(test_all main_test_all.c)
target_link_libraries(test_all PRIVATE ${LIBS_DEPEND})
