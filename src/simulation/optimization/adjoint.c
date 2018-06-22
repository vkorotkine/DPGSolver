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

#include "adjoint.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petscsys.h"

#include "macros.h"
#include "definitions_core.h"
#include "definitions_tol.h"

#include "optimization_case.h"
#include "simulation.h"

#include "multiarray.h"
#include "multiarray_constructors.h"

#include "test_case.h"
#include "solve.h"
#include "solve_implicit.h"
#include "computational_elements.h"
#include "intrusive.h"
#include "definitions_intrusive.h"

#include "volume_solver.h"


// Static function declarations ************************************************************************************* //

static Vec constructor_petsc_vec(const double*const data, const int num_vals);

static struct Multiarray_d* constructor_RHS_cmplx_stp(struct Optimization_Case *optimization_case);
static struct Multiarray_d* constructor_RHS_finite_diff(struct Optimization_Case *optimization_case);


// Interface functions ********************************************************************************************** //


void setup_adjoint(struct Optimization_Case *optimization_case){

	/*
	Setup the adjoint vectors and matrices. Compute the linearization of the objective function
	with respect to the flow (RHS) using the complex step and use the already computed linearization 
	of the residual with respect to the flow (LHS). Load all matrices and vectors into 
	PETSc structures.

	Arguments:
		optimization_case = The optimiation case data structure. The adjoint data will be 
			loaded into this structure.

	Return:
		- 
	*/

	UNUSED(constructor_RHS_finite_diff);

	// =========================
	//       RHS = [dI/dW]
	// =========================

	// Compute the RHS using the complex step
	optimization_case->RHS = constructor_RHS_cmplx_stp(optimization_case);  // destroy in solve_adjoint
	optimization_case->RHS_petsc = constructor_petsc_vec(
		(const double*const)optimization_case->RHS->data, 
		(const int)optimization_case->RHS->extents[0]
	);  // destroy in solve_adjoint


		/*
		// Testing the RHS against finite differences

		subtract_in_place_Multiarray_d(optimization_case->RHS, 
			(const struct const_Multiarray_d*)constructor_RHS_finite_diff(optimization_case));
		printf("Differernce RHS Norm: %e \n",  norm_Multiarray_d(optimization_case->RHS, "Inf"));

		exit(0);
		*/


	// =========================
	//    Chi Vector (Adjoint)
	// =========================

	// Setup the Chi vector (adjoint) to be equal to the RHS as the initial
	// vector for the iterative solver
	VecCreateSeq(MPI_COMM_WORLD,(const int)optimization_case->RHS->extents[0],&optimization_case->Chi_petsc); // keep
	VecCopy(optimization_case->RHS_petsc, optimization_case->Chi_petsc);  // destroy in solve_adjoint


	// =========================
	//      LHS = [dR/dW]^T
	// =========================

	// Obtain the LHS and load its transpose

	struct Simulation* sim = optimization_case->sim;
	struct Solver_Storage_Implicit* ssi = constructor_Solver_Storage_Implicit(sim); // destructed
	compute_rlhs_optimization(sim, ssi);

	// Copy the LHS values into a new Matrix and transpose it
	Mat LHS = ssi->A;
	Mat LHS_T;
	MatDuplicate(LHS, MAT_COPY_VALUES, &LHS_T);
	MatTranspose(LHS_T, MAT_INPLACE_MATRIX, &LHS_T);

	optimization_case->LHS_petsc = LHS_T;  // destroy in solve_adjoint


	// Free allocated structures
	destructor_Solver_Storage_Implicit(ssi);

}


void solve_adjoint(struct Optimization_Case *optimization_case){

	/*
	Solve the adjoint equation using the PETSc KSP solver. The solution will be 
	loaded into the Chi (adjoint) multiarray_d structure in optimization_case.

	Arguments:
		optimization_case = The data structure with the optimization data

	Return:
		-
	*/

	Mat LHS;
	Vec RHS;
	Vec Chi;

	LHS = optimization_case->LHS_petsc;
	RHS = optimization_case->RHS_petsc;
	Chi = optimization_case->Chi_petsc;

	/*
	printf("LHS: \n");
	MatView(LHS, PETSC_VIEWER_STDOUT_SELF);

	printf("RHS: \n");
	VecView(RHS, PETSC_VIEWER_STDOUT_SELF);

	printf("Chi: \n");
	VecView(Chi, PETSC_VIEWER_STDOUT_SELF);
	*/

	KSP ksp = NULL;
	constructor_petsc_ksp(&ksp,LHS,optimization_case->sim); // destructed

	// Get the convergence history
	//PetscReal *ksp_residual = 0;
	//PetscInt n_ksp_residual = 200;
	//PetscMalloc(n_ksp_residual*sizeof(PetscReal),&ksp_residual);
	//KSPSetResidualHistory(ksp, ksp_residual, n_ksp_residual, PETSC_FALSE);

	// Solve the system of equations
	KSPSolve(ksp, RHS, Chi);

	int iteration_ksp;
	KSPGetIterationNumber(ksp, &iteration_ksp);
	
	printf("Adjoint Solved (iter : %d ) \n", iteration_ksp);fflush(stdout);

	//printf("Get Residual History\n");
	//KSPGetResidualHistory(ksp, &ksp_residual, &n_ksp_residual);
	//int k;

	//for (k = 0; k < n_ksp_residual; k++){
	//	printf("%d %e \n", k, ksp_residual[k]);
	//}
	//exit(0);


	// Save the solution in the Chi vector (doubles)
	double *Chi_i = get_col_Multiarray_d(0, optimization_case->Chi);
	
	int *ix;
	ix = malloc((unsigned int)optimization_case->RHS->extents[0]* sizeof *ix);  // free

	for (int i = 0; i < (int)optimization_case->RHS->extents[0]; i++)
		ix[i] = i;

	VecGetValues(Chi, (int)optimization_case->RHS->extents[0], ix, Chi_i);


	// Free data structures:
	free(ix);
	destructor_Multiarray_d(optimization_case->RHS);

	// Destroy petsc structures
	KSPDestroy(&ksp);

	MatDestroy(&LHS);
	LHS = NULL;

	VecDestroy(&RHS);
	RHS = NULL;

	VecDestroy(&Chi);
	Chi = NULL;

}



// Static functions ************************************************************************************************* //

static struct Multiarray_d* constructor_RHS_cmplx_stp(struct Optimization_Case *optimization_case){

	/*
	Compute the RHS ([dI/dW], where I is the objective function) using the complex step.
	Return the results in a multiarray

	Arguments:
		optimization_case = The data structure with the optimization material

	Return:
		Multiarray_d with the RHS components
	*/

	// Preprocessing:
	// Compute the size of the RHS vector

	struct Simulation *sim_c = optimization_case->sim_c; 

	int RHS_size_ex_0;  // Number of rows (extent 0) for RHS vector

	// Use the first volume as reference to get the P value
	// NOTE: Assumes global P is constant (need to take into account adaptation
	// eventually)
	struct Solver_Volume_c* s_vol = (struct Solver_Volume_c*)sim_c->volumes->first;

	int P = s_vol->p_ref,  // P value for all elements
		nv = (int)sim_c->n_v,  // number of volumes  
		neq = NEQ_EULER,  // number of equations
		d = DIM;  // Dimension

	RHS_size_ex_0 = nv * (int)(pow(P+1, d)) * neq;

	
	// =========================
	//       RHS = [dI/dW]
	// =========================

	// Compute the RHS vector using the complex step

	struct Multiarray_d *RHS = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){RHS_size_ex_0,1});
	double *RHS_i = get_col_Multiarray_d(0, RHS);

	int RHS_index; 
	double dI_by_dW_component;

	// Loop over the complex volumes, perturb the solution vector by a complex step and compute
	// each RHS component.
	for (struct Intrusive_Link* curr = sim_c->volumes->first; curr; curr = curr->next){

		struct Solver_Volume_c *s_vol_c = (struct Solver_Volume_c*) curr;
		struct Multiarray_c* sol_coef_c = s_vol_c->sol_coef;

		// Perturb each solution vector component by the complex step and obtain the 
		// linearization terms.
		const ptrdiff_t n_col_l = compute_size(sol_coef_c->order,sol_coef_c->extents);
		for (int col_l = 0; col_l < n_col_l; col_l++) {
			sol_coef_c->data[col_l] += CX_STEP*I;
			dI_by_dW_component = cimag(optimization_case->objective_function_c(sim_c))/CX_STEP; 
			sol_coef_c->data[col_l] -= CX_STEP*I;

			// Load the value into the vector
			RHS_index = (int)s_vol_c->ind_dof+col_l;
			RHS_i[RHS_index] = dI_by_dW_component;
		}
	}

	return RHS;

}


static struct Multiarray_d* constructor_RHS_finite_diff(struct Optimization_Case *optimization_case){

	/*
	Compute the RHS ([dI/dW], where I is the objective function) using finite differences.
	Return the results in a multiarray

	Arguments:
		optimization_case = The data structure with the optimization material

	Return:
		Multiarray_d with the RHS components
	*/

	// Preprocessing:
	// Compute the size of the RHS vector

	struct Simulation *sim = optimization_case->sim; 

	int RHS_size_ex_0;  // Number of rows (extent 0) for RHS vector

	// Use the first volume as reference to get the P value
	// NOTE: Assumes global P is constant (need to take into account adaptation
	// eventually)
	struct Solver_Volume* s_vol = (struct Solver_Volume*)sim->volumes->first;

	int P = s_vol->p_ref,  // P value for all elements
		nv = (int)sim->n_v,  // number of volumes  
		neq = NEQ_EULER,  // number of equations
		d = DIM;  // Dimension

	RHS_size_ex_0 = nv * (int)(pow(P+1, d)) * neq;

	
	// =========================
	//       RHS = [dI/dW]
	// =========================

	// Compute the RHS vector using the complex step

	struct Multiarray_d *RHS = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){RHS_size_ex_0,1});
	double *RHS_i = get_col_Multiarray_d(0, RHS);

	int RHS_index; 
	double dI_by_dW_component;
	double I_0 = optimization_case->objective_function(sim);

	// Loop over the complex volumes, perturb the solution vector by a complex step and compute
	// each RHS component.
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next){

		struct Solver_Volume *s_vol = (struct Solver_Volume*) curr;
		struct Multiarray_d* sol_coef = s_vol->sol_coef;

		// Perturb each solution vector component by the complex step and obtain the 
		// linearization terms.
		const ptrdiff_t n_col_l = compute_size(sol_coef->order,sol_coef->extents);
		for (int col_l = 0; col_l < n_col_l; col_l++) {
			sol_coef->data[col_l] += FINITE_DIFF_STEP;
			dI_by_dW_component = (optimization_case->objective_function(sim) - I_0)/FINITE_DIFF_STEP;
			sol_coef->data[col_l] -= FINITE_DIFF_STEP;

			// Load the value into the vector
			RHS_index = (int)s_vol->ind_dof+col_l;
			RHS_i[RHS_index] = dI_by_dW_component;
		}
	}

	return RHS;

}


static Vec constructor_petsc_vec(const double*const data, const int num_vals){

	/*
	Construct a petsc vector from an array of data (of doubles).

	Arguments:
		data = array of doubles with the data to be loaded into the petsc vector structure.
		num_vals = The size of the data (number of elements)

	Return:
		Petsc vector with the data filled into it.
	*/

	Vec V_petsc;

	// Initialize the standard petsc vector
	VecCreateSeq(MPI_COMM_WORLD,num_vals,&V_petsc); // keep

	// ix stores the indeces of the elements in the global vector. Fill the whole vector 
	// into the Petsc vector
	PetscInt *ix;
	ix = malloc((unsigned int)num_vals* sizeof *ix);  // free
	for (int i = 0; i < num_vals; i++)
		ix[i] = i;
	
	// Set the values in the PETSc vector and assemble it
	VecSetValues(V_petsc, num_vals, ix, data, INSERT_VALUES);	
	VecAssemblyBegin(V_petsc);
	VecAssemblyEnd(V_petsc);

	free(ix);

	return V_petsc;

}


