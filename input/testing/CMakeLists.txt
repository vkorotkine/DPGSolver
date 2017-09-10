file(GLOB_RECURSE CTRL_FILES *.ctrl)
foreach(ctrl_name_full ${CTRL_FILES})
	string(REPLACE "${CMAKE_SOURCE_DIR}/input" "" ctrl_name ${ctrl_name_full})
	configure_file(${ctrl_name_full} ${PROJECT_BINARY_DIR}${ctrl_name} COPYONLY)
endforeach()

file(GLOB_RECURSE DATA_FILES *.data)
foreach(data_name_full ${DATA_FILES})
	string(REPLACE "${CMAKE_SOURCE_DIR}/input/testing/" "" data_name ${data_name_full})
	configure_file(${data_name_full} ${PROJECT_BINARY_DIR}/testing/${data_name} COPYONLY)
endforeach()

add_subdirectory(fe)
add_subdirectory(mesh)
