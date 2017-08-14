// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__mesh_connectivity_h__INCLUDED
#define DPG__mesh_connectivity_h__INCLUDED
/**	\file
 *	Add general comments for the header file here.
 */

#include "mesh_readers.h"

/// \brief Set up the mesh connectivity.
struct Mesh_Connectivity* mesh_connect
	(const struct Mesh_Data*const mesh_data ///< \ref Mesh_Data.
	);

#endif // DPG__mesh_connectivity_h__INCLUDED
