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

#include "def_templates_type_d.h"
#include "def_templates_solve_dg.h"
#include "def_templates_face_solver.h"
#include "def_templates_multiarray.h"
#include "solve_dg_T.h"
#include "undef_templates_type.h"
#include "undef_templates_solve_dg.h"
#include "undef_templates_face_solver.h"
#include "undef_templates_multiarray.h"

struct const_Matrix_d;
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
	(const struct Simulation* sim,       ///< \ref Simulation.
	 struct Solver_Storage_Implicit* ssi ///< \ref Solver_Storage_Implicit.
	);

/** \brief Set the values of \ref Solver_Storage_Implicit::row and Solver_Storage_Implicit::col based on the current
 *         volume and eq, var indices. */
void set_petsc_Mat_row_col
	(struct Solver_Storage_Implicit*const ssi, ///< \ref Solver_Storage_Implicit.
	 const struct Solver_Volume* v_l,          ///< The left \ref Solver_Volume_T.
	 const int eq,                             ///< The index of the equation.
	 const struct Solver_Volume* v_r,          ///< The right \ref Solver_Volume_T.
	 const int vr                              ///< The index of the variable.
	);

/// \brief Add lhs values to the petsc Mat at the appropriate location.
void add_to_petsc_Mat
	(const struct Solver_Storage_Implicit*const ssi, ///< \ref Solver_Storage_Implicit.
	 const struct const_Matrix_d* lhs                ///< The matrix containing the lhs data.
	);

/// \brief Version of \ref compute_flux_imbalances for the DG scheme.
void compute_flux_imbalances_dg
	(const struct Simulation*const sim ///< See brief.
	);

#endif // DPG__solve_dg_h__INCLUDED
