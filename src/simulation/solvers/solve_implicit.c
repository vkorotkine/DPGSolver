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
 */

#include "solve_implicit.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"

#include "macros.h"
#include "definitions_intrusive.h"

#include "computational_elements.h"
#include "face.h"
#include "volume_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Perform one implicit step.
 *  \return The absolute value of the maximum rhs for the current solution. */
static double implicit_step
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Display the solver progress.
static void display_progress
	(const struct Test_Case* test_case, ///< \ref Test_Case.
	 const int i_step,                  ///< The current implicit step.
	 const double max_rhs,              ///< The current maximum value of the rhs term.
	 const double max_rhs0              ///< The initial maximum value of the rhs term.
	);

/** \brief Check the exit conditions.
 *  \return `true` if exit conditions are satisfied, `false` otherwise. */
static bool check_exit
	(const struct Test_Case* test_case, ///< \ref Test_Case.
	 const double max_rhs,              ///< The current maximum value of the rhs term.
	 const double max_rhs0              ///< The initial maximum value of the rhs term.
	);

// Interface functions ********************************************************************************************** //

void solve_implicit (struct Simulation* sim)
{
	assert(sim->method == METHOD_DG); // Can be made flexible in future.

	sim->test_case->solver_method_curr = 'i';
	constructor_derived_Elements(sim,IL_ELEMENT_SOLVER_DG);       // destructed
	constructor_derived_computational_elements(sim,IL_SOLVER_DG); // destructed

	struct Test_Case* test_case = sim->test_case;

	double max_rhs0 = 0.0;
	for (int i_step = 0; ; ++i_step) {
		const double max_rhs = implicit_step(sim);
		if (i_step == 0)
			max_rhs0 = max_rhs;

		display_progress(test_case,i_step,max_rhs,max_rhs0);
		if (check_exit(test_case,max_rhs,max_rhs0))
			break;
	}

	destructor_derived_computational_elements(sim,IL_SOLVER);
	destructor_derived_Elements(sim,IL_ELEMENT);
	sim->test_case->solver_method_curr = 0;
}

// Level 0 ********************************************************************************************************** //

/** \brief Constructor for a \ref Solver_Storage_Implicit container.
 *  \return See brief. */
struct Solver_Storage_Implicit* constructor_Solver_Storage_Implicit
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Solver_Storage_Implicit container.
void destructor_Solver_Storage_Implicit
	(struct Solver_Storage_Implicit* s_store_i ///< \ref Solver_Storage_Implicit.
	);

static double implicit_step (const struct Simulation* sim)
{
	struct Solver_Storage_Implicit* s_store_i = constructor_Solver_Storage_Implicit(sim);
	const double max_rhs = compute_rlhs(sim,s_store_i);

//	Vec x; ///< Petsc Vec holding the computed solution (\f$ x = A^{-1}b \f$).
//	update_sol_coef(sim);
	destructor_Solver_Storage_Implicit(s_store_i);
EXIT_UNSUPPORTED;
	return max_rhs;
}

static void display_progress
	(const struct Test_Case* test_case, const int i_step, const double max_rhs, const double max_rhs0)
{
UNUSED(max_rhs0);
	if (!test_case->display_progress)
		return;

double emin = 1.0, emax = 0.0;
int iteration_ksp = 0;
int reason = 0;
	printf("iteration: %5d, KSP iterations (cond, reason): %5d (% .3e, %d), maxRHS: % .3e\n",
	       i_step,iteration_ksp,emax/emin,reason,max_rhs);
EXIT_ADD_SUPPORT;
}

static bool check_exit (const struct Test_Case* test_case, const double max_rhs, const double max_rhs0)
{
	bool exit_now = false;

	if (max_rhs < test_case->exit_tol_i) {
		printf("Complete: max_rhs is below the exit tolerance.\n");
		exit_now = true;
	}

	if (max_rhs/max_rhs0 < test_case->exit_ratio_i) {
		printf("Complete: max_rhs dropped by % .2e orders.\n",log10(max_rhs/max_rhs0));
		exit_now = true;
	}

	return exit_now;
}

// Level 1 ********************************************************************************************************** //

/// \brief Update \ref Solver_Volume::ind_dof and \ref Solver_Face::ind_dof.
static void update_ind_dof
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for a \ref Vector_i\* holding the 'n'umber of 'n'on-'z'ero entries in each row of the global
 *         system matrix.
 *  \return See brief. */
static struct Vector_i* constructor_nnz
	(const struct Simulation* sim ///< \ref Simulation.
	);

struct Solver_Storage_Implicit* constructor_Solver_Storage_Implicit (const struct Simulation* sim)
{
	assert(sizeof(PetscInt) == sizeof(int)); // Ensure that all is working correctly if this is removed.

	const ptrdiff_t dof = compute_dof(sim);
	update_ind_dof(sim);
	struct Vector_i* nnz = constructor_nnz(sim); // destructed

	struct Solver_Storage_Implicit* s_store_i = calloc(1,sizeof *s_store_i); // free

	MatCreateSeqAIJ(MPI_COMM_WORLD,dof,dof,0,nnz->data,&s_store_i->A); // destructed
	VecCreateSeq(MPI_COMM_WORLD,dof,&s_store_i->b);                    // destructed

	destructor_Vector_i(nnz);

	return s_store_i;
}

void destructor_Solver_Storage_Implicit (struct Solver_Storage_Implicit* s_store_i)
{
	MatDestroy(&s_store_i->A);
	VecDestroy(&s_store_i->b);

	free(s_store_i);
}

// Level 2 ********************************************************************************************************** //

/// \brief Increment the corresponding rows of `nnz` by the input number of columns.
static void increment_nnz
	(struct Vector_i* nnz,    ///< Holds the number of non-zero entries for each row.
	 const ptrdiff_t ind_dof, ///< The index of the first degree of freedom for rows to be incremented.
	 const ptrdiff_t n_row,   ///< The number of sequential rows to be incremented.
	 const ptrdiff_t n_col    ///< The increment.
	);

static void update_ind_dof (const struct Simulation* sim)
{
	ptrdiff_t dof = 0;
	switch (sim->method) {
	case METHOD_DG:
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;

			s_vol->ind_dof = dof;

			struct Multiarray_d* sol_coef = s_vol->sol_coef;
			dof += compute_size(sol_coef->order,sol_coef->extents);
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",sim->method);
		break;
	}
	assert(dof == compute_dof(sim));
}

static struct Vector_i* constructor_nnz (const struct Simulation* sim)
{
	const ptrdiff_t dof = compute_dof(sim);
	struct Vector_i* nnz = constructor_zero_Vector_i(dof); // returned

	switch (sim->method) {
	case METHOD_DG:
		// Diagonal contribution
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;

			struct Multiarray_d* sol_coef = s_vol->sol_coef;
			const ptrdiff_t size = compute_size(sol_coef->order,sol_coef->extents);
			increment_nnz(nnz,s_vol->ind_dof,size,size);
		}

		// Off-diagonal contributions
		for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
			struct Face* face = (struct Face*) curr;
			if (face->boundary)
				continue;

			struct Solver_Volume* s_vol[2] = { (struct Solver_Volume*) face->neigh_info[0].volume,
			                                   (struct Solver_Volume*) face->neigh_info[1].volume, };

			struct Multiarray_d* sol_coef[2] = { s_vol[0]->sol_coef, s_vol[1]->sol_coef, };
			const ptrdiff_t size[2] = { compute_size(sol_coef[0]->order,sol_coef[0]->extents),
			                            compute_size(sol_coef[1]->order,sol_coef[1]->extents), };

			increment_nnz(nnz,s_vol[0]->ind_dof,size[0],size[1]);
			increment_nnz(nnz,s_vol[1]->ind_dof,size[1],size[0]);
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",sim->method);
		break;
	}

	return nnz;
}

// Level 3 ********************************************************************************************************** //

static void increment_nnz (struct Vector_i* nnz, const ptrdiff_t ind_dof, const ptrdiff_t n_row, const ptrdiff_t n_col)
{
	assert(ind_dof >= 0);

	const ptrdiff_t i_max = ind_dof+n_row;
	assert(i_max <= nnz->ext_0);

	for (ptrdiff_t i = ind_dof; i < i_max; ++i)
		nnz->data[i] += n_col;
}
