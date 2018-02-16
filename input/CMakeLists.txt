file(GLOB_RECURSE GEO_FILES *.geo)
foreach(geo_name_full ${GEO_FILES})
	string(REPLACE "${CMAKE_SOURCE_DIR}/input" "" geo_name ${geo_name_full})
	configure_file(${geo_name_full} ${PROJECT_BINARY_DIR}${geo_name} COPYONLY)
endforeach()

file(GLOB_RECURSE SCRIPT_FILES *.sh)
foreach(script_name_full ${SCRIPT_FILES})
	string(REPLACE "${CMAKE_SOURCE_DIR}/input" "" script_name ${script_name_full})
	configure_file(${script_name_full} ${PROJECT_BINARY_DIR}${script_name} @ONLY)
endforeach()

add_subdirectory(meshes)
add_subdirectory(testing)
