set(SOURCE mesh.c
           mesh_readers.c
           mesh_connectivity.c
           mesh_periodic.c
           mesh_vertices.c
           mesh_geometry_cylinder_hollow_section.c
   )

add_library(Mesh ${SOURCE})
