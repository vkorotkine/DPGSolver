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

struct Matrix_d;
struct const_Matrix_d;
struct Multiarray_d;
struct Solver_Face;
struct Solver_Volume;
struct Simulation;
struct Solver_Storage_Implicit;

/** \brief Version of \ref compute_rhs for the dg method.
 *  \return See brief. */
double compute_rhs_dg
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Version of \ref compute_rlhs for the dg method.
 *  \return See brief. */
double compute_rlhs_dg
	(const struct Simulation* sim,             ///< \ref Simulation.
	 struct Solver_Storage_Implicit* s_store_i ///< \ref Solver_Storage_Implicit.
	);

/** \brief Permute the input multiarray such that its ordering is such that it is in the reference coordinates of the
 *         face cubature nodes of the opposite volume. */
void permute_Multiarray_d_fc
	(struct Multiarray_d* data,       ///< The data to be permuted.
	 const char perm_layout,          ///< Defined for \ref permute_Multiarray_d_V.
	 const int side_index_dest,       ///< The side index of the destination.
	 const struct Solver_Face* s_face ///< \ref Solver_Face.
	);

/** \brief Permute the input matrix such that its ordering is such that it is in the reference coordinates of the
 *         face cubature nodes of the opposite volume. */
void permute_Matrix_d_fc
	(struct Matrix_d* data,           ///< The data to be permuted.
	 const char perm_layout,          ///< The layout in which to permute.
	 const int side_index_dest,       ///< The side index of the destination.
	 const struct Solver_Face* s_face ///< \ref Solver_Face.
	);

/** \brief Get the pointer to the appropriate \ref DG_Solver_Element::nc_fc \ref const_Vector_i\*.
 *  \return See brief. */
const struct const_Vector_i* get_operator__nc_fc__dg
	(const int side_index_dest,       ///< Defined for \ref permute_Multiarray_d_fc.
	 const struct Solver_Face* s_face ///< Defined for \ref permute_Multiarray_d_fc.
	);

/** \brief Set the values of \ref Solver_Storage_Implicit::row and Solver_Storage_Implicit::col based on the current
 *         volume and eq, var indices. */
void set_petsc_Mat_row_col
	(struct Solver_Storage_Implicit*const s_store_i, ///< \ref Solver_Storage_Implicit.
	 const struct Solver_Volume* v_l,                ///< The left \ref Solver_Volume.
	 const int eq,                                   ///< The index of the equation.
	 const struct Solver_Volume* v_r,                ///< The right \ref Solver_Volume.
	 const int vr                                    ///< The index of the variable.
	);

/// \brief Add lhs values to the petsc Mat at the appropriate location.
void add_to_petsc_Mat
	(const struct Solver_Storage_Implicit*const s_store_i, ///< \ref Solver_Storage_Implicit.
	 const struct const_Matrix_d* lhs                      ///< The matrix containing the lhs data.
	);

#endif // DPG__solve_dg_h__INCLUDED
