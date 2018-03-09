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
#include "definitions_physics.h"
#include "definitions_test_case.h"
#include "definitions_tol.h"
#include "definitions_visualization.h"

#include "element_solver.h"
#include "face_solver.h"
#include "volume_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "computational_elements.h"
#include "intrusive.h"
#include "math_functions.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solution_euler.h"
#include "solve.h"
#include "solve_dg.h"
#include "solve_dpg.h"
#include "test_case.h"
#include "visualization.h"

// Static function declarations ************************************************************************************* //

///\{ \name Flags for whether certain outputs are enabled.
#define OUTPUT_PETSC_AB false ///< Flag for Petsc data containers.
#define OUTPUT_SOLUTION false ///< Flag for solution data on each element.
#define EXIT_ON_OUTPUT  false ///< Flag for whether the simulation should exit after outputting.
#define OUTPUT_STEP     10    ///< Iteration step at which to output the solution if enabled.
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
	(const struct Test_Case* test_case, ///< \ref Test_Case_T.
	 const double max_rhs               ///< The current maximum value of the rhs term.
	);

// Interface functions ********************************************************************************************** //

void solve_implicit (struct Simulation* sim)
{
	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	test_case->solver_method_curr = 'i';

	constructor_derived_elements_comp_elements(sim);     // destructed
	for (int i_step = 0; ; ++i_step) {
		const double max_rhs = implicit_step(i_step,sim);

		if (check_exit(test_case,max_rhs))
			break;
	}
	destructor_derived_elements_comp_elements(sim);

	test_case->solver_method_curr = 0;
}

bool check_symmetric (const struct Simulation* sim)
{
	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	const int pde_index = test_case->pde_index;

	switch (pde_index) {
	case PDE_ADVECTION:
		switch (sim->method) {
		case METHOD_DG:
			return false;
			break;
		case METHOD_DPG:
			return true;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",sim->method);
			break;
		}
		break;
	case PDE_DIFFUSION:
		return true;
		break;
	case PDE_EULER:         // fallthrough
	case PDE_NAVIER_STOKES:
		// From some previous testing, it is possible that the LHS matrix is symmetric when the test norm is
		// not a function of the solution but this will generally not be the case, notably for the quasi-optimal
		// test norm.
		return false;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",pde_index);
		break;
	}
}

void output_petsc_mat (Mat A, const char* file_name)
{
	PetscViewer viewer;

	PetscViewerASCIIOpen(PETSC_COMM_WORLD,file_name,&viewer);
	PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	MatView(A,viewer);

	PetscViewerDestroy(&viewer);
}

// Level 0 ********************************************************************************************************** //

/// \brief Output the petsc Mat/Vec to a file for visualization.
static void output_petsc_mat_vec
	(Mat A,                       ///< The petsc Mat.
	 Vec b,                       ///< The petsc Vec.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Output the solution for visualization.
static void output_solution
	(const int i_step,           ///< Defined for \ref implicit_step.
	 struct Simulation*const sim ///< Defined for \ref implicit_step.
	);

/** \brief Solve the global system of equations and update the coefficients corresponding to the dof.
 *  \return Petsc error code. */
static PetscErrorCode solve_and_update
	(const double max_rhs,                      ///< The maximum of the rhs terms.
	 const int i_step,                          ///< Defined for \ref implicit_step.
	 const struct Solver_Storage_Implicit* ssi, ///< \ref Solver_Storage_Implicit.
	 const struct Simulation* sim               ///< \ref Simulation.
	);

/** \brief Check if the pde under consideration is linear.
 *  \return `true` if yes; `false` otherwise. */
static bool check_pde_linear
	(const int pde_index ///< \ref Test_Case_T::pde_index.
	);

static void constructor_derived_elements_comp_elements (struct Simulation* sim)
{
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
}

static double implicit_step (const int i_step, const struct Simulation* sim)
{
	struct Solver_Storage_Implicit* ssi = constructor_Solver_Storage_Implicit(sim); // destructed

	printf("\tCompute rlhs.\n");
	const double max_rhs = compute_rlhs(sim,ssi);

	petsc_mat_vec_assemble(ssi);
	if (OUTPUT_PETSC_AB)
		output_petsc_mat_vec(ssi->A,ssi->b,sim);

	if (OUTPUT_SOLUTION && i_step == OUTPUT_STEP)
		output_solution(i_step,(struct Simulation*)sim);

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
		printf("Complete: max_rhs dropped by % .2e orders.\n",log10(max_rhs0/max_rhs));
		exit_now = true;
		max_rhs0 = 0.0;
	}

	if (check_pde_linear(test_case->pde_index))
		exit_now = true;

	return exit_now;
}

// Level 1 ********************************************************************************************************** //

#define N_SCHUR 4 ///< The number of sub-blocks extracted for the Schur complement.

/// \brief Container for data relating to Schur complement of the input matrix and vector.
struct Schur_Data {
	PetscInt* idx[2]; ///< The indices to be stored in the index sets.

	/** The index sets corresponding to the two regions of the square matrix to be decomposed into the Schur
	 *  complement components. Note that only two index sets are provided as the row and column sets are the same
	 *  for symmetric matrices. */
	IS is[2];

	Mat* submat;   ///< The array of sub-matrices.
	Vec subvec[2]; ///< The array of sub-vectors.
};

/// \brief Output a PETSc Vec to the file of given input name.
static void output_petsc_vec
	(Vec b,                ///< The PETSc Vec.
	 const char* file_name ///< The file name.
	);

/// \brief Output the Schur complement sub-matrices and sub-vectors to a file for visualization.
static void output_petsc_schur
	(Mat A,                       ///< The petsc Mat.
	 Vec b,                       ///< The petsc Vec.
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
 *  \return The Petsc error code. */
static PetscErrorCode constructor_petsc_ksp
	(KSP*const ksp,               ///< Pointer to the Petsc KSP.
	 Mat A,                       ///< The matrix.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Update the values of coefficients based on the computed increment.
static void update_coefs
	(Vec x,                       ///< Petsc Vec holding the solution coefficient increments.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Display the solver progress.
static void display_progress
	(const struct Test_Case* test_case, ///< \ref Test_Case_T.
	 const int i_step,                  ///< The current implicit step.
	 const double max_rhs,              ///< The current maximum value of the rhs term.
	 KSP ksp                            ///< Petsc `KSP` context.
	);

/** \brief Constructor for a \ref Schur_Data container.
 *  \return See brief. */
static struct Schur_Data* constructor_Schur_Data
	(Mat A,                       ///< The full PETSc Mat.
	 Vec b,                       ///< The full PETSc Vec.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Schur_Data container.
static void destructor_Schur_Data
	(Vec b,                        ///< The full PETSc Vec from which the sub-vectors were obtained.
	 struct Schur_Data* schur_data ///< See brief.
	);

static void output_petsc_mat_vec (Mat A, Vec b, const struct Simulation* sim)
{
	output_petsc_mat(A,"A.m");
	output_petsc_vec(b,"b.m");

	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	if (test_case->use_schur_complement)
		output_petsc_schur(A,b,sim);

	EXIT_ERROR("Disable outputting to continue");
}

static void output_solution (const int i_step, struct Simulation*const sim)
{
	UNUSED(i_step);

	output_visualization(sim,VIS_GEOM_EDGES);
	output_visualization(sim,VIS_SOLUTION);

	if (EXIT_ON_OUTPUT)
		EXIT_ERROR("Disable outputting to continue");
}

static PetscErrorCode solve_and_update
	(const double max_rhs, const int i_step, const struct Solver_Storage_Implicit* ssi, const struct Simulation* sim)
{
	KSP ksp = NULL;
	Vec x = constructor_petsc_x(ssi->b); // destructed

	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	const bool use_schur_complement = test_case->use_schur_complement;
	if (!use_schur_complement) {
		printf("\tKSP set up.\n");
		CHKERRQ(constructor_petsc_ksp(&ksp,ssi->A,sim)); // destructed
		printf("\tKSP solve.\n");
		CHKERRQ(KSPSolve(ksp,ssi->b,x));
	} else {
		printf("\tCompute Schur.\n");
		struct Schur_Data* schur_data = constructor_Schur_Data(ssi->A,ssi->b,sim); // destructed
		Mat A[2][2] = { { schur_data->submat[0], schur_data->submat[1], },
		                { schur_data->submat[2], schur_data->submat[3], } };
		Vec b[2] = { schur_data->subvec[0], schur_data->subvec[1], };

		/// \todo Run with '-info' petsc option and check that fill has the appropriate value.
		const double fill = 1.0;

		Mat A_01_11i,
		    A_01_11i_10;

		CHKERRQ(MatMatMult(A[0][1], A[1][1],MAT_INITIAL_MATRIX,fill,&A_01_11i));    // destroyed
		CHKERRQ(MatMatMult(A_01_11i,A[1][0],MAT_INITIAL_MATRIX,fill,&A_01_11i_10)); // destroyed

		CHKERRQ(MatAXPY(A[0][0],-1.0,A_01_11i_10,SAME_NONZERO_PATTERN));
		CHKERRQ(MatDestroy(&A_01_11i_10));

		CHKERRQ(VecScale(b[1],-1.0));
		CHKERRQ(MatMultAdd(A_01_11i,b[1],b[0],b[0]));
		CHKERRQ(MatDestroy(&A_01_11i));
		CHKERRQ(VecScale(b[1],-1.0));

		Vec x_sub[2];
		for (int i = 0; i < 2; ++i)
			CHKERRQ(VecGetSubVector(x,schur_data->is[i],&x_sub[i])); // restored

		printf("\tKSP set up.\n");
		CHKERRQ(constructor_petsc_ksp(&ksp,A[0][0],sim)); // destructed

		printf("\tKSP solve.\n");
		CHKERRQ(KSPSolve(ksp,b[0],x_sub[0]));

		Mat A_11i_10;
		CHKERRQ(MatMatMult(A[1][1],A[1][0],MAT_INITIAL_MATRIX,fill,&A_11i_10)); // destroyed

		CHKERRQ(MatMult(A_11i_10,x_sub[0],x_sub[1]));
		CHKERRQ(MatDestroy(&A_11i_10));
		CHKERRQ(VecScale(x_sub[1],-1.0));
		CHKERRQ(MatMultAdd(A[1][1],b[1],x_sub[1],x_sub[1]));

		for (int i = 0; i < 2; ++i)
			CHKERRQ(VecRestoreSubVector(x,schur_data->is[i],&x_sub[i]));

		destructor_Schur_Data(ssi->b,schur_data);
	}

	update_coefs(x,sim);
	destructor_petsc_x(x);

	display_progress(test_case,i_step,max_rhs,ksp);
	CHKERRQ(KSPDestroy(&ksp));
	return 0;
}

static bool check_pde_linear (const int pde_index)
{
	switch (pde_index) {
	case PDE_ADVECTION: // fallthrough
	case PDE_DIFFUSION:
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

/// \brief Update the values of \ref Solver_Volume_T::sol_coef based on the computed increment.
static void update_coef_s_v
	(Vec x,                       ///< Petsc Vec holding the solution coefficient increments.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Update the values of \ref Solver_Face_T::nf_coef based on the computed increment.
static void update_coef_nf_f
	(Vec x,                       ///< Petsc Vec holding the solution coefficient increments.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Update the values of \ref Solver_Volume_T::l_mult based on the computed increment.
static void update_coef_l_mult_v
	(Vec x,                            ///< Petsc Vec holding the solution coefficient increments.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

static void output_petsc_schur (Mat A, Vec b, const struct Simulation* sim)
{
	struct Schur_Data* schur_data = constructor_Schur_Data(A,b,sim); // destructed
	Mat* submat = schur_data->submat;

	for (int i = 0; i < N_SCHUR; ++i) {
		char mat_name[STRLEN_MIN] = { 0, };
		sprintf(mat_name,"%c%d%d%s",'A',i/2,i%2,".m");
		output_petsc_mat(submat[i],mat_name);
	}

	destructor_Schur_Data(b,schur_data);
}

static void output_petsc_vec (Vec b, const char* file_name)
{
	PetscViewer viewer;

	PetscViewerASCIIOpen(PETSC_COMM_WORLD,file_name,&viewer);
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
	VecSetFromOptions(x);
	VecSetUp(x);
	VecCopy(b,x);

	VecAssemblyBegin(x);
	VecAssemblyEnd(x);

	return x;
}

static void destructor_petsc_x (Vec x)
{
	VecDestroy(&x);
}

static PetscErrorCode constructor_petsc_ksp (KSP*const ksp, Mat A, const struct Simulation* sim)
{
	CHKERRQ(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
	CHKERRQ(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

	CHKERRQ(KSPCreate(MPI_COMM_WORLD,ksp));
	CHKERRQ(KSPSetOperators(*ksp,A,A));
	CHKERRQ(KSPSetFromOptions(*ksp));
	CHKERRQ(KSPSetComputeSingularValues(*ksp,PETSC_TRUE));

	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;

	const bool symmetric = check_symmetric(sim);
	const int solver_type_i = test_case->solver_type_i;
	switch (solver_type_i) {
	case SOLVER_I_DIRECT: {
		PC pc;
		CHKERRQ(KSPGetPC(*ksp,&pc));

		CHKERRQ(KSPSetType(*ksp,KSPPREONLY));
		CHKERRQ(KSPSetInitialGuessNonzero(*ksp,PETSC_FALSE));
		if (!symmetric)
			CHKERRQ(PCSetType(pc,PCLU));
		else
			CHKERRQ(PCSetType(pc,PCCHOLESKY));
		CHKERRQ(PCSetUp(pc));
		break;
	} case SOLVER_I_ITERATIVE:
/// \todo Potentially modify the tolerance used here based on the current residual value.
//		KSPSetTolerances(*ksp,1e-15,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
		break;
	}
	KSPSetUp(*ksp);
	return 0;
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
		update_coef_l_mult_v(x,sim);
		assert(((struct Test_Case*)sim->test_case_rc->tc)->has_2nd_order == false); // Add support.
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

	if (reason < 0) {
		printf("Warning: Petsc solver diverged with KSPConvergedReason: %d.\n",reason);
		if (reason != KSP_DIVERGED_ITS)
			EXIT_ERROR("Diverged for reason other than maximum iterations.\n");
	}

	printf("iteration: %5d, KSP iterations (cond, reason): %5d (% .3e, %d), max rhs (initial): % .3e (% .3e)\n",
	       i_step,iteration_ksp,emax/emin,reason,max_rhs,max_rhs0);
}

static struct Schur_Data* constructor_Schur_Data (Mat A, Vec b, const struct Simulation* sim)
{
	const ptrdiff_t dof[2] = { compute_dof_schur('k',sim), compute_dof_schur('v',sim), };

	struct Schur_Data* schur_data = malloc(sizeof *schur_data); // free

	PetscInt** idx = schur_data->idx;
	IS* is = schur_data->is;
	for (int i = 0; i < 2; ++i) {
		idx[i] = malloc((size_t)dof[i] * sizeof *idx[i]); // free
		const PetscInt dof_base = (PetscInt)( i == 0 ? 0 : dof[0] );
		for (int j = 0; j < dof[i]; ++j)
			idx[i][j] = dof_base+j;
		ISCreateGeneral(MPI_COMM_WORLD,(PetscInt)dof[i],idx[i],PETSC_USE_POINTER,&is[i]); // destroyed
		VecGetSubVector(b,is[i],&schur_data->subvec[i]); // restored
	}

	IS is_row[N_SCHUR] = { is[0], is[0], is[1], is[1], },
	   is_col[N_SCHUR] = { is[0], is[1], is[0], is[1], };
	MatCreateSubMatrices(A,N_SCHUR,is_row,is_col,MAT_INITIAL_MATRIX,&schur_data->submat); // destroyed

	return schur_data;
}

static void destructor_Schur_Data (Vec b, struct Schur_Data* schur_data)
{
	for (int i = 0; i < 2; ++i) {
		VecRestoreSubVector(b,schur_data->is[i],&schur_data->subvec[i]);
		ISDestroy(&schur_data->is[i]);
		free(schur_data->idx[i]);
	}
	MatDestroySubMatrices(N_SCHUR,&schur_data->submat);
	free(schur_data);
}

// Level 3 ********************************************************************************************************** //

/// \brief Update the input coefficients with the step stored in the PETSc Vec.
static void update_coef
	(const int ind_dof,                 ///< The index of the first dof associated with the coefficients.
	 struct Multiarray_d*const coef,    ///< The coefficients to be updated.
	 Vec x,                             ///< The PETSc Vec holding the updates.
	 const bool allow_relax,            ///< Flag for whether under relaxation should be allowed.
	 const struct Solver_Volume* s_vol, ///< The current \ref Solver_Volume_T.
	 const struct Simulation* sim       ///< \ref Simulation.
	);

static void update_coef_s_v (Vec x, const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol = (struct Solver_Volume*)curr;

		update_coef((int)s_vol->ind_dof,s_vol->sol_coef,x,true,s_vol,sim);
		enforce_positivity_highorder(s_vol,sim);
	}
}

static void update_coef_nf_f (Vec x, const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Solver_Face* s_face = (struct Solver_Face*)curr;

		update_coef((int)s_face->ind_dof,s_face->nf_coef,x,false,NULL,NULL);
	}
}

static void update_coef_l_mult_v (Vec x, const struct Simulation*const sim)
{
	if (!test_case_explicitly_enforces_conservation(sim))
		return;

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol = (struct Solver_Volume*)curr;

		update_coef((int)s_vol->ind_dof_constraint,s_vol->l_mult,x,false,NULL,NULL);
	}
}

// Level 4 ********************************************************************************************************** //

/** \brief Compute the under relaxation required to maintain physically correct data.
 *  \return See brief.
 *
 *  The magnitude of the under relaxation is set such that the average value of the density and pressure remain
 *  positive after the perturbation is added to the coefficients.
 */
static double compute_relaxation
	(const struct const_Multiarray_d*const coef,   ///< The coefficients.
	 const struct const_Multiarray_d*const d_coef, ///< The full perturbation to the coefficients.
	 const struct Solver_Volume*const s_vol,       ///< \ref Solver_Volume_T.
	 const struct Simulation*const sim             ///< \ref Simulation.
	);

static void update_coef
	(const int ind_dof, struct Multiarray_d*const coef, Vec x, const bool allow_relax,
	 const struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	assert(sizeof(PetscScalar) == sizeof(double)); // Use appropriate multiarray otherwise.

	const int ni = (int)compute_size(coef->order,coef->extents);
	if (ni == 0)
		return;

	PetscInt ix[ni];
	for (int i = 0; i < ni; ++i)
		ix[i] = ind_dof+i;

	PetscScalar y[ni];
	VecGetValues(x,ni,ix,y);

	double relax = 1.0;
	if (allow_relax) {
		assert(s_vol != NULL);
		const struct const_Multiarray_d* d_coef =
			constructor_move_const_Multiarray_d_d(coef->layout,coef->order,coef->extents,false,y); // destructed
		relax = compute_relaxation((struct const_Multiarray_d*)coef,(struct const_Multiarray_d*)d_coef,s_vol,sim);
		destructor_const_Multiarray_d(d_coef);
	}

	for (int i = 0; i < ni; ++i)
		coef->data[i] += relax*y[i];
}

// Level 5 ********************************************************************************************************** //

static double compute_relaxation
	(const struct const_Multiarray_d*const coef, const struct const_Multiarray_d*const d_coef,
	 const struct Solver_Volume*const s_vol, const struct Simulation*const sim)
{
	double relax = 1.0;
	if (!test_case_requires_positivity((struct Test_Case*) sim->test_case_rc->tc))
		return relax;

	const struct Volume* vol         = (struct Volume*) s_vol;
	const struct Solver_Element* s_e = (struct Solver_Element*) vol->element;

	const int p = s_vol->p_ref;
	const struct Operator* ccSB0_vs_vs = get_Multiarray_Operator(s_e->ccSB0_vs_vs,(ptrdiff_t[]){0,0,p,p});

	const char op_format = 'd';
	struct Multiarray_d* coef_b = constructor_mm_NN1_Operator_Multiarray_d(
		ccSB0_vs_vs,(struct Multiarray_d*)coef,'C',op_format,coef->order,NULL); // dest.
	const struct const_Multiarray_d* d_coef_b =
		constructor_mm_NN1_Operator_const_Multiarray_d(ccSB0_vs_vs,d_coef,'C',op_format,d_coef->order,NULL); // dest.

	const ptrdiff_t n_n  = coef_b->extents[0],
	                n_vr = coef_b->extents[1];

	ptrdiff_t* extents = (ptrdiff_t[]) { 1, n_vr, };
	double data[n_vr];
	struct Multiarray_d vr_avgs = { .layout = 'C', .order = 2, .extents = extents, .owns_data = false, .data = data, };

	do {
		add_in_place_Multiarray_d(relax,coef_b,d_coef_b);
		convert_variables(coef_b,'c','p');

		for (int vr = 0; vr < n_vr; ++vr) {
			const double*const vr_data = &coef_b->data[vr*n_n];
			vr_avgs.data[vr] = average_d(vr_data,n_n);
		}

		if (!((vr_avgs.data[0]      < EPS_PHYS) ||
		      (vr_avgs.data[n_vr-1] < EPS_PHYS)))
			break;

		convert_variables(coef_b,'p','c');
		add_in_place_Multiarray_d(-relax,coef_b,d_coef_b);
		relax /= 2.0;
	} while (relax > EPS);

	if (relax < EPS)
		EXIT_ERROR("Under relaxation approaching zero.\n");

	destructor_Multiarray_d(coef_b);
	destructor_const_Multiarray_d(d_coef_b);

	return relax;
}
