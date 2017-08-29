set	(SOURCE
	 mesh.c
	 mesh_readers.c
	 mesh_connectivity.c
	 mesh_periodic.c
	 mesh_vertices.c
	 mesh_geometry_cylinder_hollow_section.c
	)

set	(LIBS
	 General
	)

add_library(Mesh ${SOURCE})
target_link_libraries(Mesh ${LIBS})
