// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__setup_mesh_h__INCLUDED
#define DPG__setup_mesh_h__INCLUDED
/// \file

#include "Parameters.h"

/// \brief Hold data related to the mesh processing.
struct Mesh_Data {
	unsigned int n_ve,               ///< Number of vertices.
	             n_pve,              ///< Number of periodic vertices.
	             n_elem_dim[DMAX+1], ///< Number of elements of each dimension.
	             n_elem_tot,         ///< Sum of \ref n_elem_dim components

	unsigned int* pve,        ///< Holds the list of periodic vertices.
	            * elem_types, ///< Holds the list of element types (gmsh convention) as read from the mesh file.
	            * elem_tags,  /**< Holds the list of tags associated with each element.
	                           *   - 1st tag: Physical tag (associated with gmsh's physical elements); gives boundary
	                           *     condition information.
	                           *   - 2nd tag: Index of the geometry element in gmsh associated with this element; needed
	                           *     for the determination of periodic connections if present.
	                           */

	double* ve_xyz, ///< Vertex xyz locations.
};

extern void setup_mesh (void);

#endif // DPG__setup_mesh_h__INCLUDED
