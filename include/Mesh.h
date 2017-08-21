// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__set_up_mesh_h__INCLUDED
#define DPG__set_up_mesh_h__INCLUDED
/**	\file
 *	\brief Provides an interface for reading a mesh file.
 *
 *	Relevant data is returned as part of \ref Mesh_Data.
 *
 *	Supported input formats:
 *		- gmsh
 */

#include "Intrusive.h"
#include "mesh_readers.h"
#include "mesh_connectivity.h"


/// \brief Holds data output from the mesh set up.
struct Mesh {
	const struct Mesh_Data*const         mesh_data; ///< \ref Mesh_Data.
	const struct Mesh_Connectivity*const mesh_conn; ///< \ref Mesh_Connectivity.
};


/** \brief Set up the \ref Mesh.
 *
 *	To provide addtional modularity, it is possible to pass a `NULL` value for the `elements` list. This results in the
 *	list being constructed and destructed as part of the connectivity set up. However, as the base \ref Element list is
 *	used in many other modules of the code, it is generally convenient to set it up before setting up the mesh.
 *	\todo Add support for this functionality.
 */
struct Mesh* set_up_mesh
	(const char*const mesh_name_full,            ///< Defined in \ref Simulation.
	 const unsigned int d,                       ///< Defined in \ref Simulation.
	 const struct const_Intrusive_List* elements ///< The base \ref Element list.
	);

/// \brief Destructor for a \ref Mesh.
void destructor_Mesh
	(struct Mesh* mesh ///< Standard.
	);

/** \brief See return.
 *	\return The index of the first volume. */
size_t get_first_volume_index
	(const struct const_Vector_ui*const elem_per_dim, ///< Defined in \ref Conn_info.
	 const unsigned int d                             ///< Defined in \ref Conn_info.
	);


#endif // DPG__set_up_mesh_h__INCLUDED
