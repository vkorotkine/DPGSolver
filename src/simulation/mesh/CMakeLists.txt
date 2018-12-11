set	(SOURCE
	 mesh.cpp
	 mesh_readers.cpp
	 mesh_connectivity.cpp
	 mesh_periodic.cpp
	 mesh_vertices.cpp
	 mesh_geometry_cylinder_hollow_section.cpp
	)

set	(LIBS_DEPEND
	 Containers
	 General
	 Element
	)

add_library(Mesh ${SOURCE})
target_link_libraries(Mesh ${LIBS_DEPEND})
