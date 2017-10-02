find_package(MPI REQUIRED)

find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})

find_package(GSL 2.4 REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

include_directories(
	allocators
	containers
	definitions
	general
	intrusive
	simulation
	simulation/computational_elements
	simulation/element
	simulation/geometry
	simulation/mesh
	simulation/operators
	simulation/solution
	testing
	)

add_subdirectory(allocators)
add_subdirectory(containers)
add_subdirectory(general)
add_subdirectory(intrusive)
add_subdirectory(simulation)
add_subdirectory(testing)
