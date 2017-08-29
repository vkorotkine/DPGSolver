include_directories(
	allocators
	constants
	containers
	general
	intrusive
	simulation
	simulation/mesh
	testing
	)

add_subdirectory(allocators)
add_subdirectory(containers)
add_subdirectory(general)
add_subdirectory(intrusive)
add_subdirectory(simulation)
add_subdirectory(testing)
