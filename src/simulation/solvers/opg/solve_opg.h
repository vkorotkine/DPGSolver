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

#ifndef DPG__solve_opg_h__INCLUDED
#define DPG__solve_opg_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used to solve for the solution using the 'o'ptimal 'p'etrov
 *         'g'alerkin method.
 */

#include "def_templates_type_d.h"
#include "def_templates_solve_opg.h"
#include "solve_opg_T.h"
#include "undef_templates_type.h"
#include "undef_templates_solve_opg.h"

struct OPG_Solver_Volume;

/** \brief Version of \ref compute_rlhs for the opg method.
 *  \return See brief. */
double compute_rlhs_opg
	(const struct Simulation*const sim,       ///< Standard.
	 struct Solver_Storage_Implicit*const ssi ///< Standard.
	);

/** \brief Set the values of \ref Solver_Storage_Implicit::row and Solver_Storage_Implicit::col based on the input
 *         volumes and eq, var indices. */
void set_petsc_Mat_row_col_opg
	(struct Solver_Storage_Implicit*const ssi, ///< \ref Solver_Storage_Implicit.
	 const struct OPG_Solver_Volume*const v_l, ///< The left \ref Solver_Volume_T.
	 const int eq,                             ///< The index of the equation.
	 const struct OPG_Solver_Volume*const v_r, ///< The right \ref Solver_Volume_T.
	 const int vr                              ///< The index of the variable.
	);

#endif // DPG__solve_opg_h__INCLUDED
