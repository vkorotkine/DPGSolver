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

#ifndef DPG__test_complex_compute_all_rhs_dpg_h__INCLUDED
#define DPG__test_complex_compute_all_rhs_dpg_h__INCLUDED
/** \file
 *  \brief Provides `complex` versions of functions defined in \ref compute_all_rlhs_dpg.h.
 */

struct Simulation;
struct Solver_Storage_Implicit;
struct Solver_Volume;
struct Complex_DPG_Solver_Volume;

/// \brief Version of \ref compute_all_rlhs_dpg computing only complex rhs terms.
void compute_all_rhs_dpg_c
	(const struct Complex_DPG_Solver_Volume* c_dpg_s_vol, ///< The \ref Complex_DPG_Solver_Volume.
	 struct Solver_Storage_Implicit* ssi,                 ///< See brief.
	 const struct Simulation* sim                         ///< See brief.
	);

/** \brief See \ref constructor_sol_vc_interp.
 *  \return Standard. */
const struct const_Multiarray_c* constructor_sol_vc_dpg_c
	(const struct Solver_Volume* s_vol, ///< See brief.
	 const struct Simulation* sim       ///< See brief.
	);

#endif // DPG__test_complex_compute_all_rhs_dpg_h__INCLUDED
