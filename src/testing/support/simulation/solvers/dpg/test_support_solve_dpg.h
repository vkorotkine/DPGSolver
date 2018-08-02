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

#ifndef DPG__test_support_solve_dpg_h__INCLUDED
#define DPG__test_support_solve_dpg_h__INCLUDED
/** \file
 *  \brief Provides supporting functions for testing of the dpg solver.
 */

struct Simulation;
struct Solver_Storage_Implicit;

/// \brief Compute the lhs matrix using the complex step method for the DG scheme.
void compute_lhs_cmplx_step_dpg
	(const struct Simulation* sim,       ///< \ref Simulation.
	 struct Solver_Storage_Implicit* ssi ///< \ref Solver_Storage_Implicit.
	);

/** \brief Version of \ref constructor_Gen_Eig_Data_inf_sup for the dpg method.
 *  \return See brief. */
const struct Gen_Eig_Data* constructor_Gen_Eig_Data_inf_sup_dpg
	(struct Simulation*const sim ///< Standard.
	 );

#endif // DPG__test_support_solve_dpg_h__INCLUDED
