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

#ifndef DPG__solve_dg_h__INCLUDED
#define DPG__solve_dg_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used to solve for the solution using the 'd'iscontinuous 'g'alerkin
 *         method.
 */

struct Multiarray_d;
struct Solver_Face;
struct Simulation;

/** \brief Version of \ref compute_rhs for the dg method.
 *  \return See brief. */
double compute_rhs_dg
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Permute the input multiarray such that its ordering is such that it is in the reference coordinates of the
 *         face cubature nodes of the opposite volume. */
void permute_Multiarray_d_fc
	(struct Multiarray_d* data,       ///< The data to be permuted.
	 const char perm_layout,          ///< Defined for \ref permute_Multiarray_d_V.
	 const int side_index_dest,       ///< The side index of the destination.
	 const struct Solver_Face* s_face ///< \ref Solver_Face.
	);

#endif // DPG__solve_dg_h__INCLUDED
