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
#include "face_solver.h"
#include "volume_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "solve_dg.h"
#include "solve_dpg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

///\{ \name Flag for whether the petsc data containers should be output to a file.
#define OUTPUT_PETSC_AB false
///\}

/// \brief Constructor for the derived element and computational element lists.
static void constructor_derived_elements_comp_elements
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for the derived element and computational element lists.
static void destructor_derived_elements_comp_elements
	(struct Simulation* sim ///< \ref Simulation.
	);

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

/** \brief Constructor for a \ref Vector_i\* holding the 'n'umber of 'n'on-'z'ero entries in each row of the global
 *         system matrix.
 *  \return See brief. */
static struct Vector_i* constructor_nnz
	(const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void solve_implicit (struct Simulation* sim)
{
	sim->test_case->solver_method_curr = 'i';

	constructor_derived_elements_comp_elements(sim); // destructed
	for (int i_step = 0; ; ++i_step) {
		const double max_rhs = implicit_step(i_step,sim);

		if (check_exit(sim->test_case,max_rhs))
			break;
	}
	destructor_derived_elements_comp_elements(sim);

	sim->test_case->solver_method_curr = 0;
}

struct Solver_Storage_Implicit* constructor_Solver_Storage_Implicit (const struct Simulation* sim)
{
	assert(sizeof(PetscInt) == sizeof(int)); // Ensure that all is working correctly if this is removed.

	const ptrdiff_t dof = compute_dof(sim);
	update_ind_dof(sim);
	struct Vector_i* nnz = constructor_nnz(sim); // destructed

	struct Solver_Storage_Implicit* ssi = calloc(1,sizeof *ssi); // free

	MatCreateSeqAIJ(MPI_COMM_WORLD,dof,dof,0,nnz->data,&ssi->A); // destructed
	VecCreateSeq(MPI_COMM_WORLD,dof,&ssi->b);                    // destructed

	destructor_Vector_i(nnz);

	return ssi;
}

void destructor_Solver_Storage_Implicit (struct Solver_Storage_Implicit* ssi)
{
	MatDestroy(&ssi->A);
	VecDestroy(&ssi->b);

	free(ssi);
}

void petsc_mat_vec_assemble (struct Solver_Storage_Implicit* ssi)
{
	MatAssemblyBegin(ssi->A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ssi->A,MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(ssi->b);
	VecAssemblyEnd(ssi->b);
}

void increment_nnz (struct Vector_i* nnz, const ptrdiff_t ind_dof, const ptrdiff_t n_row, const ptrdiff_t n_col)
{
	assert(ind_dof >= 0);

	const ptrdiff_t i_max = ind_dof+n_row;
	assert(i_max <= nnz->ext_0);

	for (ptrdiff_t i = ind_dof; i < i_max; ++i)
		nnz->data[i] += n_col;
}

bool check_symmetric (const struct Simulation* sim)
{
	switch (sim->method) {
	case METHOD_DG: {
		const int pde_index = sim->test_case->pde_index;
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
		break;
	} case METHOD_DPG:
		return true;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
}

// Level 0 ********************************************************************************************************** //

/// \brief Output the petsc Mat/Vec to a file for visualization.
static void output_petsc_mat_vec
	(Mat A,                       ///< The petsc Mat.
	 Vec b,                       ///< The petsc Vec.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Solve the global system of equations and update the coefficients corresponding to the dof.
static void solve_and_update
	(const double max_rhs,                      ///< The maximum of the rhs terms.
	 const int i_step,                          ///< Defined for \ref implicit_step.
	 const struct Solver_Storage_Implicit* ssi, ///< \ref Solver_Storage_Implicit.
	 const struct Simulation* sim               ///< \ref Simulation.
	);

/** \brief Check if the pde under consideration is linear.
 *  \return `true` if yes; `false` otherwise. */
static bool check_pde_linear
	(const int pde_index ///< \ref Test_Case::pde_index.
	);

static void constructor_derived_elements_comp_elements (struct Simulation* sim)
{
	constructor_derived_Elements(sim,IL_ELEMENT_SOLVER); // destructed
	switch (sim->method) {
	case METHOD_DG:
		constructor_derived_Elements(sim,IL_ELEMENT_SOLVER_DG);       // destructed
		constructor_derived_computational_elements(sim,IL_SOLVER_DG); // destructed
		break;
	case METHOD_DPG:
		constructor_derived_Elements(sim,IL_ELEMENT_SOLVER_DPG);       // destructed
		constructor_derived_computational_elements(sim,IL_SOLVER_DPG); // destructed
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
}

static void destructor_derived_elements_comp_elements (struct Simulation* sim)
{
	destructor_derived_computational_elements(sim,IL_SOLVER);
	destructor_derived_Elements(sim,IL_ELEMENT_SOLVER);
	destructor_derived_Elements(sim,IL_ELEMENT);
}

static double implicit_step (const int i_step, const struct Simulation* sim)
{
	struct Solver_Storage_Implicit* ssi = constructor_Solver_Storage_Implicit(sim); // destructed

	const double max_rhs = compute_rlhs(sim,ssi);

	petsc_mat_vec_assemble(ssi);
	if (OUTPUT_PETSC_AB)
		output_petsc_mat_vec(ssi->A,ssi->b,sim);

	solve_and_update(max_rhs,i_step,ssi,sim);
	destructor_Solver_Storage_Implicit(ssi);

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

static struct Vector_i* constructor_nnz (const struct Simulation* sim)
{
	struct Vector_i* nnz = NULL;
	switch (sim->method) {
	case METHOD_DG:  nnz = constructor_nnz_dg(sim);  break;
	case METHOD_DPG: nnz = constructor_nnz_dpg(sim); break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",sim->method);
		break;
	}
	return nnz;
}

// Level 1 ********************************************************************************************************** //

#define N_SCHUR 4 ///< The number of sub-blocks extracted for the Schur complement.

/// \brief Output the Schur complement petsc Mat sub-matrices to a file for visualization.
static void output_petsc_mat_schur
	(Mat A,                       ///< The petsc Mat.
	 const struct Simulation* sim ///< \ref Simulation.
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

/// \brief Update the values of coefficients based on the computed increment.
static void update_coefs
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

/** \brief Constructor for the Schur complement sub-matrices of the input matrix.
 *  \return See brief. */
static Mat* constructor_schur_sub_matrices
	(Mat A,                       ///< The input PETSc Mat.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for the Schur complement sub-matrices.
static void destructor_schur_sub_matrices
	(Mat* submat ///< The sub-matrices.
	);

static void output_petsc_mat_vec (Mat A, Vec b, const struct Simulation* sim)
{
	PetscViewer viewer;

	PetscViewerASCIIOpen(PETSC_COMM_WORLD,"A.m",&viewer);
	PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	MatView(A,viewer);

	PetscViewerASCIIOpen(PETSC_COMM_WORLD,"b.m",&viewer);
	PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	VecView(b,viewer);

	PetscViewerDestroy(&viewer);

	if (sim->test_case->use_schur_complement)
		output_petsc_mat_schur(A,sim);

	EXIT_ERROR("Disable outputting to continue");
}

static void solve_and_update
	(const double max_rhs, const int i_step, const struct Solver_Storage_Implicit* ssi, const struct Simulation* sim)
{
	KSP ksp = NULL;
	Vec x   = NULL;

	const bool use_schur_complement = sim->test_case->use_schur_complement;
	if (!use_schur_complement) {
		x   = constructor_petsc_x(ssi->b);       // destructed
		ksp = constructor_petsc_ksp(ssi->A,sim); // destructed
		KSPSolve(ksp,ssi->b,x);
	} else {
		Mat* submat = constructor_schur_sub_matrices(ssi->A,sim); // destructed

		PetscInt sub_nnz[N_SCHUR];
		for (int i = 0; i < N_SCHUR; ++i) {
			MatInfo matinfo;
			MatGetInfo(submat[i],MAT_GLOBAL_MAX,&matinfo);

			sub_nnz[i] = matinfo.nz_used;
printf("%d\n",sub_nnz[i]);
		}
		const PetscReal fill = sub_nnz[3]/((PetscReal)(sub_nnz[0]+sub_nnz[1]+sub_nnz[2]));

// fill - expected fill as ratio of nnz(D)/(nnz(A) + nnz(B)+nnz(C))
		Mat D;
		MatMatMatMult(submat[1],submat[3],submat[2],MAT_INITIAL_MATRIX,fill,&D); // destroyed

		MatDestroy(&D);

		destructor_schur_sub_matrices(submat);

		EXIT_ADD_SUPPORT;
	}

	update_coefs(x,sim);
	destructor_petsc_x(x);

	display_progress(sim->test_case,i_step,max_rhs,ksp);
	destructor_petsc_ksp(ksp);
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

/// \brief Update the values of \ref Solver_Volume::sol_coef based on the computed increment.
static void update_coef_s_v
	(Vec x,                       ///< Petsc Vec holding the solution coefficient increments.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Update the values of \ref Solver_Face::nf_coef based on the computed increment.
static void update_coef_nf_f
	(Vec x,                       ///< Petsc Vec holding the solution coefficient increments.
	 const struct Simulation* sim ///< \ref Simulation.
	);

static void output_petsc_mat_schur (Mat A, const struct Simulation* sim)
{
	Mat* submat = constructor_schur_sub_matrices(A,sim); // destructed

	PetscViewer viewer;
	for (int i = 0; i < N_SCHUR; ++i) {
		char mat_name[STRLEN_MIN] = { 0, };
		sprintf(mat_name,"%c%d%d%s",'A',i/2,i%2,".m");

		PetscViewerASCIIOpen(PETSC_COMM_WORLD,mat_name,&viewer);
		PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
		MatView(submat[i],viewer);
	}
	PetscViewerDestroy(&viewer);

	destructor_schur_sub_matrices(submat);
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

	const bool symmetric = check_symmetric(sim);
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

static void update_coefs (Vec x, const struct Simulation* sim)
{
	switch (sim->method) {
	case METHOD_DG:
		update_coef_s_v(x,sim);
		break;
	case METHOD_DPG:
		update_coef_s_v(x,sim);
		update_coef_nf_f(x,sim);
		assert(sim->test_case->has_2nd_order == false); // Add support.
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
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

static Mat* constructor_schur_sub_matrices (Mat A, const struct Simulation* sim)
{
	const ptrdiff_t dof_fv[2] = { compute_dof_schur('f',sim), compute_dof_schur('v',sim), };

	PetscInt idx_f[dof_fv[0]];
	for (int i = 0; i < dof_fv[0]; ++i)
		idx_f[i] = i;

	PetscInt idx_v[dof_fv[1]];
	for (int i = 0; i < dof_fv[1]; ++i)
		idx_v[i] = dof_fv[0]+i;

	IS is_fv[2];
	ISCreateGeneral(MPI_COMM_WORLD,dof_fv[0],idx_f,PETSC_USE_POINTER,&is_fv[0]);
	ISCreateGeneral(MPI_COMM_WORLD,dof_fv[1],idx_v,PETSC_USE_POINTER,&is_fv[1]);

	IS is_row[N_SCHUR] = { is_fv[0], is_fv[0], is_fv[1], is_fv[1], },
	   is_col[N_SCHUR] = { is_fv[0], is_fv[1], is_fv[0], is_fv[1], };

	Mat* submat;
	MatCreateSubMatrices(A,N_SCHUR,is_row,is_col,MAT_INITIAL_MATRIX,&submat); // returned
	ISDestroy(&is_fv[0]);
	ISDestroy(&is_fv[1]);

	return submat;
}

static void destructor_schur_sub_matrices (Mat* submat)
{
	MatDestroySubMatrices(N_SCHUR,&submat);
}

// Level 3 ********************************************************************************************************** //

/// \brief Update the input coefficients with the step stored in the PETSc Vec.
static void update_coef
	(const int ind_dof,              ///< The index of the first dof associated with the coefficients.
	 struct Multiarray_d*const coef, ///< The coefficients to be updated.
	 Vec x                           ///< The PETSc Vec holding the updates.
	);

static void update_coef_s_v (Vec x, const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol = (struct Solver_Volume*)curr;

		update_coef(s_vol->ind_dof,s_vol->sol_coef,x);
//		enforce_positivity_highorder(coef);
	}
}

static void update_coef_nf_f (Vec x, const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Solver_Face* s_face = (struct Solver_Face*)curr;

		update_coef(s_face->ind_dof,s_face->nf_coef,x);
	}
}

// Level 4 ********************************************************************************************************** //

static void update_coef (const int ind_dof, struct Multiarray_d*const coef, Vec x)
{
		const int ni = compute_size(coef->order,coef->extents);

		PetscInt ix[ni];
		for (int i = 0; i < ni; ++i)
			ix[i] = ind_dof+i;

		PetscScalar y[ni];
		VecGetValues(x,ni,ix,y);

		for (int i = 0; i < ni; ++i)
			coef->data[i] += y[i];
}
