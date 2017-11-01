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
#include "definitions_test_case.h"

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

///\{ \name Flag for whether the petsc data containers should be output to a file.
#define OUTPUT_PETSC_AB false
///\}

/** \brief Perform one implicit step.
 *  \return The absolute value of the maximum rhs for the current solution. */
static double implicit_step
	(const int i_step,            ///< The current implicit step.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Check the exit conditions.
 *  \return `true` if exit conditions are satisfied, `false` otherwise. */
static bool check_exit
	(const struct Test_Case* test_case, ///< \ref Test_Case.
	 const double max_rhs               ///< The current maximum value of the rhs term.
	);

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

// Interface functions ********************************************************************************************** //

void solve_implicit (struct Simulation* sim)
{
	assert(sim->method == METHOD_DG); // Can be made flexible in future.

	sim->test_case->solver_method_curr = 'i';
	constructor_derived_Elements(sim,IL_ELEMENT_SOLVER_DG);       // destructed
	constructor_derived_computational_elements(sim,IL_SOLVER_DG); // destructed

	struct Test_Case* test_case = sim->test_case;

	for (int i_step = 0; ; ++i_step) {
		const double max_rhs = implicit_step(i_step,sim);

		if (check_exit(test_case,max_rhs))
			break;
	}

	destructor_derived_computational_elements(sim,IL_SOLVER);
	destructor_derived_Elements(sim,IL_ELEMENT);
	sim->test_case->solver_method_curr = 0;
}

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

void petsc_mat_vec_assemble (struct Solver_Storage_Implicit* s_store_i)
{
	MatAssemblyBegin(s_store_i->A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(s_store_i->A,MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(s_store_i->b);
	VecAssemblyEnd(s_store_i->b);
}

// Level 0 ********************************************************************************************************** //

/// \brief Output the petsc Mat/Vec to a file for visualization.
static void output_petsc_mat_vec
	(Mat A, ///< The petsc Mat.
	 Vec b  ///< The petsc Vec.
	);

/** \brief Constructor for the `x` petsc Vec in "Ax = b" which is initialized to `b`.
 *  \return See brief. */
static Vec constructor_petsc_x
	(Vec b ///< \ref Solver_Storage_Implicit::b.
	);

/// \brief Destructor for the `x` petsc Vec.
static void destructor_petsc_x
	(Vec x ///< Standard.
	);

/** \brief Constructor for a petsc `KSP` context.
 *  \return See brief. */
static KSP constructor_petsc_ksp
	(Mat A,                       ///< The matrix.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a petsc `KSP` context.
static void destructor_petsc_ksp
	(KSP ksp ///< Standard.
	);

/// \brief Update the values of \ref Solver_Volume::sol_coef based on the computed increment.
static void update_sol_coef
	(Vec x,                       ///< Petsc Vec holding the solution coefficient increments.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Display the solver progress.
static void display_progress
	(const struct Test_Case* test_case, ///< \ref Test_Case.
	 const int i_step,                  ///< The current implicit step.
	 const double max_rhs,              ///< The current maximum value of the rhs term.
	 KSP ksp                            ///< Petsc `KSP` context.
	);

/// \brief Increment the corresponding rows of `nnz` by the input number of columns.
static void increment_nnz
	(struct Vector_i* nnz,    ///< Holds the number of non-zero entries for each row.
	 const ptrdiff_t ind_dof, ///< The index of the first degree of freedom for rows to be incremented.
	 const ptrdiff_t n_row,   ///< The number of sequential rows to be incremented.
	 const ptrdiff_t n_col    ///< The increment.
	);

/** \brief Check if the pde under consideration is linear.
 *  \return `true` if yes; `false` otherwise. */
static bool check_pde_linear
	(const int pde_index ///< \ref Test_Case::pde_index.
	);

static double implicit_step (const int i_step, const struct Simulation* sim)
{
	struct Solver_Storage_Implicit* s_store_i = constructor_Solver_Storage_Implicit(sim);

	const double max_rhs = compute_rlhs(sim,s_store_i);

	petsc_mat_vec_assemble(s_store_i);
	if (OUTPUT_PETSC_AB)
		output_petsc_mat_vec(s_store_i->A,s_store_i->b);

	Vec x = constructor_petsc_x(s_store_i->b);

	KSP ksp = constructor_petsc_ksp(s_store_i->A,sim);
	KSPSolve(ksp,s_store_i->b,x);
	destructor_Solver_Storage_Implicit(s_store_i);

	update_sol_coef(x,sim);
	destructor_petsc_x(x);

	display_progress(sim->test_case,i_step,max_rhs,ksp);
	destructor_petsc_ksp(ksp);

	return max_rhs;
}

static bool check_exit (const struct Test_Case* test_case, const double max_rhs)
{
	bool exit_now = false;

	if (max_rhs < test_case->exit_tol_i) {
		printf("Complete: max_rhs is below the exit tolerance.\n");
		exit_now = true;
	}

	static double max_rhs0 = 0.0;
	if (max_rhs0 == 0.0)
		max_rhs0 = max_rhs;

	if (max_rhs/max_rhs0 < test_case->exit_ratio_i) {
		printf("Complete: max_rhs dropped by % .2e orders.\n",log10(max_rhs/max_rhs0));
		exit_now = true;
	}

	if (check_pde_linear(test_case->pde_index))
		exit_now = true;

	return exit_now;
}

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

// Level 1 ********************************************************************************************************** //

static bool check_symmetric
	(const int pde_index ///< \ref Test_Case::pde_index
	);

static void output_petsc_mat_vec (Mat A, Vec b)
{
	PetscViewer viewer;

	PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat_output.m",&viewer);
	PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	MatView(A,viewer);

	PetscViewerASCIIOpen(PETSC_COMM_WORLD,"vec_output.m",&viewer);
	PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	VecView(b,viewer);

	PetscViewerDestroy(&viewer);
}

static Vec constructor_petsc_x (Vec b)
{
	PetscInt dof = 0;
	VecGetSize(b,&dof);

	Vec x = NULL;
	VecCreateSeq(MPI_COMM_WORLD,dof,&x); // destructed
	VecCopy(b,x);

	VecAssemblyBegin(x);
	VecAssemblyEnd(x);

	return x;
}

static void destructor_petsc_x (Vec x)
{
	VecDestroy(&x);
}

static KSP constructor_petsc_ksp (Mat A, const struct Simulation* sim)
{
	KSP ksp;
	KSPCreate(MPI_COMM_WORLD,&ksp);

	PC pc;

	KSPSetOperators(ksp,A,A);
	KSPSetComputeSingularValues(ksp,PETSC_TRUE);
	KSPGetPC(ksp,&pc);

	const bool symmetric = check_symmetric(sim->test_case->pde_index);
	const int solver_type_i = sim->test_case->solver_type_i;
	switch (solver_type_i) {
	case SOLVER_I_DIRECT:
		KSPSetType(ksp,KSPPREONLY);
		if (!symmetric)
			PCSetType(pc,PCLU);
		else
			PCSetType(pc,PCCHOLESKY);
		break;
	case SOLVER_I_ITER_DEF:
/// \todo Modify the tolerance used here.
		KSPSetTolerances(ksp,1e-15,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
		KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

		if (!symmetric) {
#if 1
			KSPSetType(ksp,KSPGMRES);
			KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization);
			KSPGMRESSetRestart(ksp,60); // Default: 30
			PCSetType(pc,PCILU);
#else
			KSPSetType(ksp,KSPRICHARDSON);
			PCSetType(pc,PCSOR);
			PCSORSetSymmetric(pc,SOR_SYMMETRIC_SWEEP);
#endif
		} else {
			KSPSetType(ksp,KSPCG);
			PCSetType(pc,PCILU);
		}
		PCFactorSetLevels(pc,1); // Cannot use MatOrdering with 0 fill
		PCFactorSetMatOrderingType(pc,MATORDERINGRCM);
		break;
	}
	KSPSetUp(ksp);
	PCSetUp(pc);

	return ksp;
}

static void destructor_petsc_ksp (KSP ksp)
{
	KSPDestroy(&ksp);
}


static void update_sol_coef (Vec x, const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol = (struct Solver_Volume*)curr;
		const int ind_dof             = s_vol->ind_dof;
		struct Multiarray_d* sol_coef = s_vol->sol_coef;

		const int ni = compute_size(sol_coef->order,sol_coef->extents);

		PetscInt ix[ni];
		for (int i = 0; i < ni; ++i)
			ix[i] = ind_dof+i;

		PetscScalar y[ni];
		VecGetValues(x,ni,ix,y);

		const double alpha = 1.0;
//		const double alpha = compute_underRelax(s_vol,y,sim);
//		enforce_positivity_highorder(sol_coef);

		for (int i = 0; i < ni; ++i)
			sol_coef->data[i] += alpha*y[i];
	}
}

static void display_progress (const struct Test_Case* test_case, const int i_step, const double max_rhs, KSP ksp)
{
	if (!test_case->display_progress)
		return;

	static double max_rhs0 = 0.0;
	if (i_step == 0)
		max_rhs0 = max_rhs;

	KSPConvergedReason reason;
	KSPGetConvergedReason(ksp,&reason);

	PetscInt iteration_ksp;
	KSPGetIterationNumber(ksp,&iteration_ksp);

	PetscReal emax = 0.0, emin = 0.0;
	KSPComputeExtremeSingularValues(ksp,&emax,&emin);

	printf("iteration: %5d, KSP iterations (cond, reason): %5d (% .3e, %d), max rhs (initial): % .3e (% .3e)\n",
	       i_step,iteration_ksp,emax/emin,reason,max_rhs,max_rhs0);
}

static bool check_pde_linear (const int pde_index)
{
	switch (pde_index) {
	case PDE_ADVECTION: // fallthrough
	case PDE_POISSON:
		return true;
		break;
	case PDE_EULER: // fallthrough
	case PDE_NAVIER_STOKES:
		return false;
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",pde_index);
		break;
	}
}

// Level 2 ********************************************************************************************************** //

static bool check_symmetric (const int pde_index)
{
	switch (pde_index) {
	case PDE_POISSON:
		return true;
		break;
	case PDE_ADVECTION:     // fallthrough
	case PDE_EULER:         // fallthrough
	case PDE_NAVIER_STOKES:
		return false;
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",pde_index);
		break;
	}
}

static void increment_nnz (struct Vector_i* nnz, const ptrdiff_t ind_dof, const ptrdiff_t n_row, const ptrdiff_t n_col)
{
	assert(ind_dof >= 0);

	const ptrdiff_t i_max = ind_dof+n_row;
	assert(i_max <= nnz->ext_0);

	for (ptrdiff_t i = ind_dof; i < i_max; ++i)
		nnz->data[i] += n_col;
}
