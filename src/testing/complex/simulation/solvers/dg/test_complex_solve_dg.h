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

#ifndef DPG__test_complex_solve_dg_h__INCLUDED
#define DPG__test_complex_solve_dg_h__INCLUDED
/** \file
 *  \brief Provides `complex` versions of functions defined in \ref solve_dg.h.
 */

struct Multiarray_c;
struct Simulation;
struct Solver_Storage_Implicit;
struct Solver_Face;

/// \brief Perturb the initial solution for the DG method.
void perturb_solution_dg
	(const struct Simulation* sim ///< Defined for \ref perturb_solution_fptr.
	);

/// \brief Set the initial solution for the complex computational elements from the real solution for the DG scheme.
void set_initial_solution_complex_dg
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Compute the lhs matrix using the complex step method for the DG scheme.
void compute_lhs_cmplx_step_dg
	(const struct Simulation* sim,             ///< \ref Simulation.
	 struct Solver_Storage_Implicit* s_store_i ///< \ref Solver_Storage_Implicit.
	);

/// \brief `complex` version of \ref permute_Multiarray_d_fc.
void permute_Multiarray_c_fc
	(struct Multiarray_c* data,       ///< See brief.
	 const char perm_layout,          ///< See brief.
	 const int side_index_dest,       ///< See brief.
	 const struct Solver_Face* s_face ///< See brief.
	);

#endif // DPG__test_complex_solve_dg_h__INCLUDED
