find_package(PythonLibs 3 REQUIRED)
find_package(NumPy REQUIRED)

if(PYTHONLIBS_FOUND)
	add_custom_target(meshes COMMAND make -f ${PROJECT_BINARY_DIR}/meshes/Makefile_meshes meshes)

	configure_file(${PROJECT_SOURCE_DIR}/input/meshes/Makefile_meshes.in
	               ${PROJECT_BINARY_DIR}/meshes/Makefile_meshes @ONLY)
endif()
