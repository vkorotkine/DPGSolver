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
/*	unsigned int n_ve,               ///< Number of vertices.
	             n_pve,              ///< Number of periodic vertices.
	             n_elem_dim[DMAX+1], ///< Number of elements of each dimension.
	             n_elem_tot,         ///< Sum of \ref n_elem_dim components
*/

//	unsigned int* pve,        ///< Holds the list of periodic vertices.
//	            * elem_types, ///< Holds the list of element types (gmsh convention) as read from the mesh file.
//	            * elem_tags,  /**< Holds the list of tags associated with each element.
/*	                           *   - 1st tag: Physical tag (associated with gmsh's physical elements); gives boundary
	                           *     condition information.
	                           *   - 2nd tag: Index of the geometry element in gmsh associated with this element; needed
	                           *     for the determination of periodic connections if present.
	                           */

//	double* ve_xyz, ///< Vertex xyz locations.
};

/// \brief Holds data relating to the mesh connectivity.
struct Mesh_Connectivity {
};


///< \brief Set up the \ref Mesh.
struct Mesh* set_up_mesh
	(const char*const mesh_name_full, ///< Defined in \ref Simulation.
	 const unsigned int d             ///< Defined in \ref Simulation.
	);


#endif // DPG__set_up_mesh_h__INCLUDED
