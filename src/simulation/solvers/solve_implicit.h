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
struct Vector_i;

/// \brief Solve for the solution using an implicit solver.
void solve_implicit
	(struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for a \ref Solver_Storage_Implicit container.
 *  \return See brief. */
struct Solver_Storage_Implicit* constructor_Solver_Storage_Implicit
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Solver_Storage_Implicit container.
void destructor_Solver_Storage_Implicit
	(struct Solver_Storage_Implicit* s_store_i ///< \ref Solver_Storage_Implicit.
	);

/// \brief Assemble \ref Solver_Storage_Implicit::A and \ref Solver_Storage_Implicit::b.
void petsc_mat_vec_assemble
	(struct Solver_Storage_Implicit* s_store_i ///< \ref Solver_Storage_Implicit.
	);

/// \brief Increment the corresponding rows of `nnz` by the input number of columns.
void increment_nnz
	(struct Vector_i* nnz,    ///< Holds the number of non-zero entries for each row.
	 const ptrdiff_t ind_dof, ///< The index of the first degree of freedom for rows to be incremented.
	 const ptrdiff_t n_row,   ///< The number of sequential rows to be incremented.
	 const ptrdiff_t n_col    ///< The increment.
	);

/** \brief Check whether the matrix under consideration is symmetric based on \ref Simulation::method and
 *         \ref Test_Case::pde_index.
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
