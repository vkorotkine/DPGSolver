add_subdirectory(containers)
add_subdirectory(element)
add_subdirectory(general)
add_subdirectory(simulation)

set	(SOURCE
	 test_support.cpp
	 test_support_fe_init.cpp
	)

set	(LIBS_DEPEND
	)

add_library(Test_Support STATIC ${SOURCE})
target_link_libraries(Test_Support ${LIBS_DEPEND})
