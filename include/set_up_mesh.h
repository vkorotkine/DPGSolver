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

#include "mesh_readers.h"


/// \brief Holds data output from the mesh set up.
struct Mesh {
	const struct Mesh_Data*const         mesh_data; ///< \ref Mesh_Data.
	const struct Mesh_Connectivity*const mesh_conn; ///< \ref Mesh_Connectivity.
};

/// \brief Holds data relating to the mesh connectivity.
struct Mesh_Connectivity {
};


///< \brief Set up the \ref Mesh.
struct Mesh* set_up_mesh
	(const char*const mesh_name_full, ///< Defined in \ref Simulation.
	 const unsigned int d             ///< Defined in \ref Simulation.
	);

///< \brief Destructor for a \ref Mesh.
void destructor_Mesh
	(struct Mesh* mesh ///< Standard.
	);


#endif // DPG__set_up_mesh_h__INCLUDED
