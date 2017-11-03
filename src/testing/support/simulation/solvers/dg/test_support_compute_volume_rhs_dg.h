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

#ifndef DPG__test_support_compute_volume_rhs_dg_h__INCLUDED
#define DPG__test_support_compute_volume_rhs_dg_h__INCLUDED
/** \file
 *  \brief Provides support functions for testing relating to computation of rhs terms for the 'd'iscontinuous
 *         'g'alerkin solver method.
 */

struct Simulation;
struct Intrusive_List;

/// \brief Version of \ref compute_volume_rlhs_dg computing only complex rhs terms.
void compute_volume_rhs_dg_c
	(const struct Simulation* sim,  ///< \ref Simulation.
	 struct Intrusive_List* volumes ///< The list of volumes for which to compute the rhs term.
	);

#endif // DPG__test_support_compute_volume_rhs_dg_h__INCLUDED
