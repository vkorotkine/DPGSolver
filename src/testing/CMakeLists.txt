include_directories(
	integration
	)

add_subdirectory(integration)

# Executable generation
add_executable(test_all main_test.c)
target_link_libraries(test_all Test_Integration)
