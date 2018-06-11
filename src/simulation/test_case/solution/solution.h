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

#ifndef DPG__solution_h__INCLUDED
#define DPG__solution_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used for solution specification (initialization).
 */

#include "def_templates_type_d.h"
#include "solution_T.h"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "solution_T.h"
#include "undef_templates_type.h"

/** \brief Identical to \ref constructor_sol_v_T but for \ref Solver_Volume_T::rhs.
 *  \return See brief. */
struct Multiarray_d* constructor_rhs_v
	(const struct Simulation* sim, ///< See brief.
	 struct Solver_Volume* s_vol,  ///< See brief.
	 const char node_kind          ///< See brief.
	);

/** \brief Constructor for the solution coefficients in the Bezier basis.
 *  \return See brief. */
struct Multiarray_d* constructor_s_coef_bezier
	(const struct Solver_Volume*const s_vol, ///< \ref Solver_Volume_T.
	 const struct Simulation*const sim       ///< \ref Simulation.
	);

#endif // DPG__solution_h__INCLUDED
