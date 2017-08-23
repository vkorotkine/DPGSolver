// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__mesh_vertices_h__INCLUDED
#define DPG__mesh_vertices_h__INCLUDED
/**	\file
 *	\todo Add general comments for the header file here.
 */

#include "Mesh.h"

/// \brief Container for additional information relating to the mesh vertices.
struct Mesh_Vertices {
	const struct const_Vector_i*const ve_curved;   ///< Flags for which vertices are located on curved boundaries.
	const struct const_Vector_i*const ve_boundary; ///< Flags for which vertices are located on boundaries.
	const struct const_Vector_i*const ve_bc;       /**< Holds one of the boundary conditions associated with the
	                                                *   vertices. The value is ambiguous (but not important) for
	                                                *   vertices touching two separate boundaries. */
};

/** \brief Obtain additional information about the mesh vertices and correct their position (if necessary).
 *
 *	The correction of the initial vertex position is required because of the limited precision of the output of the mesh
 *	file. This is only performed for curved domains which are not mapped.
 */
struct Mesh_Vertices* mesh_process_vertices
	(const struct Mesh*const mesh,                ///< \ref Mesh.
	 const struct const_Intrusive_List* elements, ///< The base \ref Element list.
	 const struct Mesh_Input*const mesh_input     ///< \ref Mesh_Input.
	);

/// \brief Destructor for \ref Mesh_Vertices.
void destructor_Mesh_Vertices
	(struct Mesh_Vertices* mesh_vert ///< Standard.
	);

#endif // DPG__mesh_vertices_h__INCLUDED
