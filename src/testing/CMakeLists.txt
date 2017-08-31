include_directories(
	integration
	)

add_subdirectory(
	integration
	)

# Library generation
set	(SOURCE
	 test_base.c
	)

add_library(Test_Base ${SOURCE})

# Executable generation
add_executable(test_all main_test.c)
target_link_libraries(test_all Test_Base Test_Integration)
