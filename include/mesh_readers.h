// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__mesh_readers_h__INCLUDED
#define DPG__mesh_readers_h__INCLUDED
/**	\file
 *	Provides interface to mesh readers.
 */

#include <stddef.h>

///\{ Definitions for the gmsh physical numbering convention.
#define BC_PERIODIC_MIN 50
#define PERIODIC_XL     51
#define PERIODIC_XR     52
#define PERIODIC_YL     53
#define PERIODIC_YR     54
#define PERIODIC_ZL     55
#define PERIODIC_ZR     56
///\}

/// \brief Holds data read from the mesh file.
struct Mesh_Data {
	const struct const_Matrix_d*const nodes; ///< The xyz coordinates of the mesh elements.

	const struct const_Vector_ui*const       elem_types;      ///< The list of element types.
	const struct const_Matrix_ui*const       elem_tags;       ///< The list of element tags.
	const struct const_Multiarray_Vector_ui*const node_nums; ///< The list of node numbers for the elements.

	const struct const_Matrix_ui*const periodic_corr; ///< The periodic entity correspondence.
};

/// \brief Set the mesh data from the input mesh file.
struct Mesh_Data* mesh_reader
	(const char*const mesh_name_full, ///< Mesh file name.
	 const unsigned int d             ///< Dimension.
	);

#endif // DPG__mesh_readers_h__INCLUDED
