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

#include "def_templates_type_dc.h"
#include "def_templates_solve_dg.h"
#include "def_templates_face_solver.h"
#include "def_templates_multiarray.h"
#include "solve_dg_T.h"
#include "undef_templates_type.h"
#include "undef_templates_solve_dg.h"
#include "undef_templates_face_solver.h"
#include "undef_templates_multiarray.h"

struct Simulation;
struct Solver_Storage_Implicit;

/// \brief Perturb the initial solution for the DG method.
void perturb_solution_dg
	(const struct Simulation* sim ///< Defined for \ref perturb_solution_fptr.
	);

/// \brief Compute the lhs matrix using the complex step method for the DG scheme.
void compute_lhs_cmplx_step_dg
	(const struct Simulation* sim,       ///< \ref Simulation.
	 struct Solver_Storage_Implicit* ssi ///< \ref Solver_Storage_Implicit.
	);

#endif // DPG__test_complex_solve_dg_h__INCLUDED
