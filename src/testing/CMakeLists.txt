# Library generation
set	(SOURCE
	 test_base.c
	)

add_library(Test_Base ${SOURCE})


include_directories(
	unit
	integration
	)

add_subdirectory(unit)
add_subdirectory(integration)

# Executable generation
add_executable(test_all main_test_all.c)
target_link_libraries(test_all Test_Base Test_Unit Test_Integration)
