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
/** \file
 *  \brief Provides the interface to templated functions used to solve for the solution using the
 *         'd'iscontinuous 'g'alerkin method.
 */

#include "def_templates_solve_dg.h"
#include "def_templates_face_solver.h"
#include "def_templates_multiarray.h"

struct Multiarray_T;
struct Solver_Face_T;
struct Simulation;

/// \brief Version of \ref update_ind_dof_T for the dg method.
void update_ind_dof_dg_T
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Version of \ref constructor_nnz for the dg method.
 *  \return See brief. */
struct Vector_i* constructor_nnz_dg_T
	(const struct Simulation* sim ///< \ref Simulation.
	);

#include "undef_templates_solve_dg.h"
#include "undef_templates_face_solver.h"
#include "undef_templates_multiarray.h"
