// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__mesh_readers_h__INCLUDED
#define DPG__mesh_readers_h__INCLUDED
/**	\file
 *	Provides the interface to mesh reader containers and functions.
 */

#include <stddef.h>

#include "Multiarray.h"
#include "Matrix.h"
#include "Vector.h"

/**\{ Definitions for the gmsh physical numbering convention.
 *	\todo Move this to where it is needed.
 */
#define BC_PERIODIC_MIN 50
#define PERIODIC_XL     51
#define PERIODIC_XR     52
#define PERIODIC_YL     53
#define PERIODIC_YR     54
#define PERIODIC_ZL     55
#define PERIODIC_ZR     56
///\}


/**\{ \name Definitions for the gmsh geometry numbering convention.
 *	\todo Move this to where it is needed.
 */
#define GMSH_XLINE_MIN  1001
#define GMSH_YLINE_MIN  2001
#define GMSH_ZLINE_MIN  3001
#define GMSH_XYFACE_MIN 4001
#define GMSH_XZFACE_MIN 5001
#define GMSH_YZFACE_MIN 6001
#define GMSH_XYZVOL_MIN 7001
///\}


/** \brief Holds data read from the mesh file.
 *
 *	A detailed description of the variables constructed from a Gmsh .msh file can be found in the
 *	[File formats][gmsh_ff] section of the Gmsh manual.
 *
 *	<!-- References: -->
 *	[gmsh_ff]: http://gmsh.info/doc/texinfo/gmsh.html#File-formats
 */
struct Mesh_Data {
	const struct const_Matrix_d*const nodes; ///< The xyz coordinates of the mesh elements.

	const struct const_Vector_ui*const            elem_types; ///< The list of element types.
	const struct const_Matrix_ui*const            elem_tags;  /**< The list of element tags.
	                                                           *   The 1st tag gives boundary condition information.
	                                                           *   The 2nd tag gives periodic connectivity information.
	                                                           */
	const struct const_Multiarray_Vector_ui*const node_nums;  ///< The list of node numbers for the elements.

	const struct const_Matrix_ui*const periodic_corr; ///< The periodic entity correspondence.
};

/// \brief Destructor for \ref Mesh_Data\*.
void destructor_Mesh_Data
	(struct Mesh_Data* mesh_data ///< Standard.
	);

/// \brief Set the mesh data from the input mesh file.
struct Mesh_Data* mesh_reader
	(const char*const mesh_name_full, ///< Mesh file name.
	 const unsigned int d             ///< Dimension.
	);

#endif // DPG__mesh_readers_h__INCLUDED
