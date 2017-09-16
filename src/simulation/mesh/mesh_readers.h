/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */

#ifndef DPG__mesh_readers_h__INCLUDED
#define DPG__mesh_readers_h__INCLUDED
/**	\file
 *	\brief Provides the interface to mesh reader containers and functions.
 */

#include <stddef.h>

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

	const struct const_Vector_i*const elem_types;           ///< The list of element types.
	const struct const_Matrix_i*const elem_tags;            /**< The list of element tags.
	                                                          *  The 1st tag gives boundary condition information.
	                                                          *  The 2nd tag gives periodic connectivity information.
	                                                          */
	const struct const_Multiarray_Vector_i*const node_nums; ///< The list of node numbers for the elements.

	/** The periodic entity correspondence.
	 *  Currently used only as an indicator for whether periodic boundaries are present.
	 *  Potentially change to simple boolean flag after testing. */
	const struct const_Matrix_i*const periodic_corr;
};

/** \brief Constructor for the \ref Mesh_Data from the input mesh file.
 *	\return Standard. */
struct Mesh_Data* constructor_Mesh_Data
	(const char*const mesh_name_full, ///< Mesh file name.
	 const int d                      ///< Dimension.
	);

/// \brief Destructor for \ref Mesh_Data\*.
void destructor_Mesh_Data
	(struct Mesh_Data* mesh_data ///< Standard.
	);

#endif // DPG__mesh_readers_h__INCLUDED
