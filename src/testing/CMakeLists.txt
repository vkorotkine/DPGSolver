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
	complex/containers
	complex/simulation/operators
	complex/simulation/solvers/dg
	complex/simulation/solvers/numerical_flux
	complex/simulation/test_case
	complex/simulation/test_case/boundary
	complex/simulation/test_case/flux
	integration
	support
	support/containers
	support/element
	support/general
	support/intrusive
	support/simulation/computational_elements
	unit
	)

add_subdirectory(complex)
add_subdirectory(integration)
add_subdirectory(support)
add_subdirectory(unit)


# Executable generation
set	(LIBS_DEPEND
	 Test_Base
	 Test_Support
	 Test_Unit 
	 Test_Integration
	)

add_executable(test_all main_test_all.c)
target_link_libraries(test_all PRIVATE ${LIBS_DEPEND})
