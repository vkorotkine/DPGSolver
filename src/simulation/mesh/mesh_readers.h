// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__mesh_readers_h__INCLUDED
#define DPG__mesh_readers_h__INCLUDED
/**	\file
 *	Provides the interface to mesh reader containers and functions.
 */

#include <stddef.h>

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

/** \brief Holds data read from the mesh file.
 *
 *	A detailed description of the variables constructed from a Gmsh .msh file can be found in the
 *	[File formats][gmsh_ff] section of the Gmsh manual.
 *
 *	<!-- References: -->
 *	[gmsh_ff]: http://gmsh.info/doc/texinfo/gmsh.html#File-formats
 */
struct Mesh_Data {
	const int d;           ///< The dimension.
	const ptrdiff_t ind_v; ///< The index of the first volume in the list of elements.

	const struct const_Vector_i*const elem_per_dim; ///< The number of elements per dimension.

	const struct const_Matrix_d*const nodes; ///< The xyz coordinates of the mesh elements (the mesh vertices).

	const struct const_Vector_i*const elem_types;            ///< The list of element types.
	const struct const_Matrix_i*const elem_tags;             /**< The list of element tags.
	                                                           *   The 1st tag gives boundary condition information.
	                                                           *   The 2nd tag gives periodic connectivity information.
	                                                           */
	const struct const_Multiarray_Vector_i*const node_nums;  ///< The list of node numbers for the elements.

	/** The periodic entity correspondence.
	 * Currently used only as an indicator for whether periodic boundaries are present.
	 * \todo Change to simple boolean flag after testing. */
	const struct const_Matrix_i*const periodic_corr;
};

/// \brief Destructor for \ref Mesh_Data\*.
void destructor_Mesh_Data
	(struct Mesh_Data* mesh_data ///< Standard.
	);

/// \brief Set the mesh data from the input mesh file.
struct Mesh_Data* mesh_reader
	(const char*const mesh_name_full, ///< Mesh file name.
	 const int d                      ///< Dimension.
	);

#endif // DPG__mesh_readers_h__INCLUDED
