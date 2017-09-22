// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__mesh_periodic_h__INCLUDED
#define DPG__mesh_periodic_h__INCLUDED
/**	\file
 *	Provides the interface to mesh periodic containers and functions.
 */

#include "mesh_connectivity.h"

/** \brief Correct the face vertex correspondence if the mesh is periodic.
 *
 *	The correction modifies the vertices of the "slave" faces to those of the "master" faces. The correspondence is
 *	established by checking the equivalence of the mesh node coordinates in the appropriate directions.
 */
void correct_f_ve_for_periodic
	(const struct Mesh_Data*const mesh_data, ///< The \ref Mesh_Data.
	 struct Conn_info*const conn_info        ///< The \ref Conn_info.
	);

#endif // DPG__mesh_periodic_h__INCLUDED
