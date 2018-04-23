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

#ifndef DPG__approximate_nearest_neighbor_h__INCLUDED
#define DPG__approximate_nearest_neighbor_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions relating to approximate nearest neighbor computation.
 *
 *  \note The only reason for the reimplementation of this functionality was a lack of c libraries with the desired
 *        features; it was determined to be simpler to do the implementation than to maintain a c wrapper for an
 *        existing c++ library. If this project is ever converted to c++, please consider directly using an available
 *        library which would almost certainly be more efficient and more extensively tested.
 */

struct const_Vector_i;
struct const_Matrix_d;

/** \brief Constructor for the \ref const_Vector_T\* (`int`) of 'a'pproximate 'n'earest 'n'eighbor indices for the input
 *         background and search node lists.
 *  \return See brief. */
const struct const_Vector_i* constructor_ann_indices
	(const struct const_Matrix_d*const nodes_b, ///< Background nodes.
	 const struct const_Matrix_d*const nodes_s  ///< Search nodes.
	);

#endif // DPG__approximate_nearest_neighbor_h__INCLUDED
