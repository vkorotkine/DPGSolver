file(GLOB_RECURSE CTRL_FILES *.ctrl)
foreach(ctrl_name_full ${CTRL_FILES})
	string(REPLACE "${CMAKE_SOURCE_DIR}/input" "" ctrl_name ${ctrl_name_full})
	configure_file(${ctrl_name_full} ${PROJECT_BINARY_DIR}${ctrl_name} COPYONLY)
endforeach()

add_subdirectory(mesh)
add_subdirectory(fe)
