# Library generation
set	(SOURCE
	 test_base.c
	)

add_library(Test_Base ${SOURCE})


include_directories(
	support
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
