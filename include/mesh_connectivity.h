// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__mesh_connectivity_h__INCLUDED
#define DPG__mesh_connectivity_h__INCLUDED
/**	\file
 *	Provides the interface to mesh connectivity containers and functions.
 */

#include "mesh_readers.h"
#include "Intrusive.h"

/// \brief Holds data relating to the mesh connectivity.
struct Mesh_Connectivity {
};

/// \brief Set up the mesh connectivity.
struct Mesh_Connectivity* mesh_connect
	(const struct Mesh_Data*const mesh_data, ///< \ref Mesh_Data.
	 const struct Intrusive_List* elements   ///< The base \ref Element list.
	);

#endif // DPG__mesh_connectivity_h__INCLUDED
