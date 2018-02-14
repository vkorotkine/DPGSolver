/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */

#ifndef DPG__mesh_h__INCLUDED
#define DPG__mesh_h__INCLUDED
/** \file
 *  \brief Provides the interface for reading a mesh file and 'converting to'/'computing' useful information.
 *
 *  Relevant data is returned as part of \ref Mesh.
 *
 *  Supported input formats:
 *  - gmsh (linear elements only).
 *
 *  The \ref Mesh_Vertices container is used to assign the `boundary` and `curved` flags to the \ref Volume and \ref
 *  Face elements in the domain. The vertex information, while generally not used for the solver, can be used to set up
 *  information relating to mesh edges (faces and volumes, by association).
 */

#include <stdbool.h>
#include <stddef.h>

struct Mesh_Input;
struct Mesh;
struct const_Vector_i;
struct const_Intrusive_List;

/// \brief Container for the required input information for the mesh processing.
struct Mesh_Input {
	const int d,           ///< The dimension.
	          domain_type; /**< The type of domain. Vertex position correction is performed for
	                        *   `domain_type` == \ref DOM_BLENDED. */

	/** Flag for whether the mesh vertices should be unrealistically corrected to lie on the input domain boundary to
	 * within a very small tolerance. See \ref mesh_vertices.h for additional discussion of this issue. */
	const bool mesh_unrealistic;

	const char* mesh_name_full; ///< Name of the mesh file (including the full path and file extension).
	const char* geom_name;      ///< Name of the base geometry to be used for the domain.
	const char* geom_spec;      ///< Additional specifications for the geometry.
};

/// \brief Container for the data output from the mesh set up.
struct Mesh {
	const struct Mesh_Data*const         mesh_data; ///< \ref Mesh_Data.
	const struct Mesh_Connectivity*const mesh_conn; ///< \ref Mesh_Connectivity.
	const struct Mesh_Vertices*const     mesh_vert; ///< \ref Mesh_Vertices.
};

/** \brief Constructor for a \ref Mesh.
 *  \return Standard.
 *
 *  To provide addtional modularity, it is possible to pass a `NULL` value for the `elements` list. This results in the
 *  list being constructed and destructed as part of the connectivity set up. However, as the base \ref Element list is
 *  used in many other modules of the code, it is generally convenient to set it up before setting up the mesh.
 */
struct Mesh* constructor_Mesh
	(const struct Mesh_Input* mesh_input,        ///< \ref Mesh_Input.
	 const struct const_Intrusive_List* elements ///< The base \ref Element list.
	);

/// \brief Destructor for a \ref Mesh.
void destructor_Mesh
	(struct Mesh* mesh ///< Standard.
	);

/** \brief See return.
 *  \return The index of the first volume. */
ptrdiff_t get_first_volume_index
	(const struct const_Vector_i*const elem_per_dim, ///< Defined in \ref Conn_info.
	 const int d                                     ///< Defined in \ref Conn_info.
	);

#endif // DPG__mesh_h__INCLUDED
