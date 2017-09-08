find_package(MPI REQUIRED)
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})

include_directories(
	allocators
	constants
	containers
	general
	intrusive
	simulation
	simulation/mesh
	simulation/geometry
	simulation/solution
	testing
	)

add_subdirectory(allocators)
add_subdirectory(containers)
add_subdirectory(general)
add_subdirectory(intrusive)
add_subdirectory(simulation)
add_subdirectory(testing)
