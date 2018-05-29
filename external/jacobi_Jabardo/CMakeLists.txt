find_package(GSL 2.4 REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

add_subdirectory(jacobi-0.9.2)
#target_include_directories(Jacobi_GSL PUBLIC jacobi-0.9.2/src)
