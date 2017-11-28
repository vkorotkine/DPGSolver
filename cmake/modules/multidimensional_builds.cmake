foreach(dim RANGE 1 3)
	string(REGEX REPLACE "_([1-3]+)D" "_${dim}D" BIN_PATH_${dim}D "${CMAKE_BINARY_DIR}/bin/")
endforeach()
