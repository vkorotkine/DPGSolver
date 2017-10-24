find_package(MPI REQUIRED)
include_directories(${MPI_C_INCLUDE_PATH})

find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})

find_package(GSL 2.4 REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

# Currently including all directories here so that `#include`s need not specify relative paths.
include_directories(
	allocators
	containers
	definitions
	element
	general
	intrusive
	processing/post
	simulation
	simulation/computational_elements
	simulation/geometry
	simulation/mesh
	simulation/operators
	simulation/solvers
	simulation/solvers/dg
	simulation/solvers/numerical_flux
	simulation/test_case
	simulation/test_case/compute_error
	simulation/test_case/flux
	simulation/test_case/solution
	simulation/test_case/solution/periodic_vortex
	simulation/test_case/solution/supersonic_vortex
	testing
	)

add_subdirectory(allocators)
add_subdirectory(containers)
add_subdirectory(element)
add_subdirectory(general)
add_subdirectory(intrusive)
add_subdirectory(processing)
add_subdirectory(simulation)
add_subdirectory(testing)
