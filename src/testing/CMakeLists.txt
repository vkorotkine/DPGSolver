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
	complex/containers
	complex/general
	complex/simulation/computational_elements
	complex/simulation/geometry
	complex/simulation/operators
	complex/simulation/solvers
	complex/simulation/solvers/dg
	complex/simulation/solvers/dpg
	complex/simulation/solvers/numerical_flux
	complex/simulation/test_case
	complex/simulation/test_case/boundary
	complex/simulation/test_case/flux
	complex/simulation/test_case/solution
	complex/simulation/test_case/solution/advection
	complex/simulation/test_case/solution/diffusion
	complex/simulation/test_case/solution/euler
	complex/simulation/test_case/solution/navier_stokes
	complex/templates
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
