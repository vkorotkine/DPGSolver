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

#ifndef DPG__inverse_mapping_h__INCLUDED
#define DPG__inverse_mapping_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions relating to inverse geometry mapping functionality (i.e. given physical
 *         coordinates, compute the associated reference coordinates).
 */

struct Matrix_d;
struct const_Matrix_d;

/** \brief Constructor for a \ref const_Matrix_T\* holding the values of the reference element coordinates corresponding
           to the input physical coordinates.
 *  \return See brief. */
struct Matrix_d* constructor_inverse_mapping_mutable
	(const int e_type,                         ///< \ref Element::type.
	 const struct const_Matrix_d*const xyz_ve, ///< The vertex xyz coordinates.
	 const struct const_Matrix_d*const xyz     ///< The xyz coordinates.
	);

/** \brief `const` version of \ref constructor_inverse_mapping_mutable.
 *  \return See brief. */
const struct const_Matrix_d* constructor_inverse_mapping
	(const int e_type,                         ///< See brief.
	 const struct const_Matrix_d*const xyz_ve, ///< See brief.
	 const struct const_Matrix_d*const xyz     ///< See brief.
	);

/** \brief Constructor for the matrix of basis functions corresponding to the p1 standard reference element vertices
 *         evaluated at the input reference coordinates.
 *  \return See brief. */
const struct const_Matrix_d* constructor_basis_std_p1
	(const int e_type,                         ///< The element type.
	 const struct const_Matrix_d*const rst_std ///< The standard reference element coordinates.
	);

#endif // DPG__inverse_mapping_h__INCLUDED
