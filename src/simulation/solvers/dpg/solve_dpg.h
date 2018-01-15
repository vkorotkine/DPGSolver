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

#ifndef DPG__solve_dpg_h__INCLUDED
#define DPG__solve_dpg_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used to solve for the solution using the 'd'iscontinuous 'p'etrov
 *         'g'alerkin method.
 */

#include "def_templates_type_d.h"
#include "def_templates_solve_dpg.h"
#include "solve_dpg_T.h"
#include "undef_templates_type.h"
#include "undef_templates_solve_dpg.h"

struct Simulation;
struct Solver_Storage_Implicit;

/** \brief Version of \ref compute_rlhs for the dpg method.
 *  \return See brief. */
double compute_rlhs_dpg
	(const struct Simulation* sim,       ///< See brief.
	 struct Solver_Storage_Implicit* ssi ///< See brief.
	);

/// \brief Version of \ref compute_flux_imbalances for the DPG scheme.
void compute_flux_imbalances_dpg
	(struct Simulation*const sim ///< \ref Simulation.
	);

#endif // DPG__solve_dpg_h__INCLUDED
