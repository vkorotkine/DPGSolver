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

#ifndef DPG__solve_implicit_h__INCLUDED
#define DPG__solve_implicit_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used to solve for the solution using implicit procedures.
 */

#include <stddef.h>
#include <stdbool.h>
#include "petscmat.h"

struct Simulation;
struct Solver_Storage_Implicit;
struct Vector_i;

/// \brief Solve for the solution using an implicit solver.
void solve_implicit
	(struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Check whether the matrix under consideration is symmetric based on \ref Simulation::method and
 *         \ref Test_Case_T::pde_index.
 *  \return `true` if symmetric; `false` otherwise. */
bool check_symmetric
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Output a PETSc Mat to the file of given input name.
void output_petsc_mat
	(Mat A,                ///< The PETSc Mat.
	 const char* file_name ///< The file name.
	);

#endif // DPG__solve_implicit_h__INCLUDED
