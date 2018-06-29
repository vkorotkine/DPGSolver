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

#include "simulation.h"

#include "multiarray.h"
#include "matrix.h"
#include "matrix_constructors.h"
#include "matrix_math.h"

#include "test_case.h"
#include "solve.h"
#include "solve_implicit.h"
#include "computational_elements.h"
#include "intrusive.h"
#include "definitions_intrusive.h"
#include "volume_solver.h"

#include "optimization_case.h"
#include "objective_functions.h"


// Static function declarations ************************************************************************************* //


static Vec constructor_petsc_vec(const double*const data, const int num_vals);


/** \brief Compute the RHS ([dI/dW], where I is the functional) using the complex step.
 * Store the results into the the adjoint_data structure in the RHS Matrix_d data structure.
 */
static void set_RHS_cmplx_stp(
	struct Adjoint_Data *adjoint_data, ///< The adjoint_data data structure to load data into
	struct Simulation *sim_c, ///< The complex simulation object to be used for the complex step
	functional_fptr_c functional_c ///< The functional function pointer (complex version for complex step)
	);


//static struct Multiarray_d* constructor_RHS_finite_diff(struct Optimization_Case *optimization_case);


// Interface functions ********************************************************************************************** //

struct Adjoint_Data* constructor_Adjoint_Data(struct Optimization_Case *optimization_case){

	struct Adjoint_Data* adjoint_data = calloc(1,sizeof *adjoint_data); // returned


	struct Simulation *sim = optimization_case->sim; 

	int num_res_eqs = 0;  // The total number of residual equations (and flow variables)
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;
		num_res_eqs += (int)compute_size(s_vol->sol_coef->order,s_vol->sol_coef->extents);
	}

	// Allocate empty RHS and Chi (adjoint) matrices
	adjoint_data->RHS = constructor_empty_Matrix_d('C', num_res_eqs, 1);
	adjoint_data->Chi = constructor_empty_Matrix_d('C', num_res_eqs, 1);

	return adjoint_data;
}


void destructor_Adjoint_Data (struct Adjoint_Data* adjoint_data){

	destructor_Matrix_d(adjoint_data->RHS);
	destructor_Matrix_d(adjoint_data->Chi);

	free((void*)adjoint_data);

}


void setup_adjoint(struct Adjoint_Data* adjoint_data, struct Simulation *sim, struct Simulation *sim_c,
	functional_fptr functional, functional_fptr_c functional_c){

	//UNUSED(constructor_RHS_finite_diff); 	// MSB: TODO: Possibly remove this
	UNUSED(sim);
	UNUSED(functional);


 	// Set the RHS terms using the complex step
	set_RHS_cmplx_stp(adjoint_data, sim_c, functional_c); 

	// Load the RHS into the Chi matrix (initial)
	for(int i = 0; i < (int)adjoint_data->RHS->ext_0; i++)
		adjoint_data->Chi->data[i] = adjoint_data->RHS->data[i];


// // Testing the RHS against finite differences
// subtract_in_place_Multiarray_d(optimization_case->RHS, 
// 	(const struct const_Multiarray_d*)constructor_RHS_finite_diff(optimization_case));
// printf("Differernce RHS Norm: %e \n",  norm_Multiarray_d(optimization_case->RHS, "Inf"));
// exit(0);


}


void solve_adjoint(struct Adjoint_Data* adjoint_data, struct Simulation *sim){

	Mat LHS_petsc;
	Vec RHS_petsc;
	Vec Chi_petsc;

	// ==================================
	//     Setup PETSc Data Structures
	// ==================================

	// RHS and Chi:
	RHS_petsc = constructor_petsc_vec(	(const double*const)adjoint_data->RHS->data, 
										(const int)adjoint_data->RHS->ext_0);  // free
	Chi_petsc = constructor_petsc_vec(	(const double*const)adjoint_data->Chi->data, 
										(const int)adjoint_data->Chi->ext_0); // free

	// LHS:
	struct Solver_Storage_Implicit* ssi = constructor_Solver_Storage_Implicit(sim); // free
	compute_rlhs_optimization(sim, ssi);

	// Copy the LHS values into a new Matrix and transpose it
	LHS_petsc = ssi->A;
	Mat LHS_T_petsc;  // free
	MatDuplicate(LHS_petsc, MAT_COPY_VALUES, &LHS_T_petsc);
	MatTranspose(LHS_T_petsc, MAT_INPLACE_MATRIX, &LHS_T_petsc);


	// ==================================
	//        Perform KSP Solve
	// ==================================

	// Setup the KSP Solve and then perform the solve
	KSP ksp = NULL;
	constructor_petsc_ksp(&ksp,LHS_T_petsc,sim); // destructed

	// Get the convergence history
	//PetscReal *ksp_residual = 0;
	//PetscInt n_ksp_residual = 200;
	//PetscMalloc(n_ksp_residual*sizeof(PetscReal),&ksp_residual);
	//KSPSetResidualHistory(ksp, ksp_residual, n_ksp_residual, PETSC_FALSE);

	// Solve the system of equations
	KSPSolve(ksp, RHS_petsc, Chi_petsc);

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


	// ==================================
	//          Save Solution
	// ==================================

	// Save the solution in the Chi vector (doubles)
	double *Chi_i = get_col_Matrix_d(0, adjoint_data->Chi);
	
	int *ix;
	ix = malloc((unsigned int)adjoint_data->Chi->ext_0* sizeof *ix);  // free
	for (int i = 0; i < (int)adjoint_data->Chi->ext_0; i++)
		ix[i] = i;

	VecGetValues(Chi_petsc, (int)adjoint_data->Chi->ext_0, ix, Chi_i);



	// Free data structures:
	free(ix);
	destructor_Solver_Storage_Implicit(ssi);

	// Destroy petsc structures
	KSPDestroy(&ksp);

	MatDestroy(&LHS_T_petsc);
	LHS_T_petsc = NULL;

	VecDestroy(&RHS_petsc);
	RHS_petsc = NULL;

	VecDestroy(&Chi_petsc);
	Chi_petsc = NULL;

}



// Static functions ************************************************************************************************* //

static void set_RHS_cmplx_stp(struct Adjoint_Data *adjoint_data, struct Simulation *sim_c, functional_fptr_c functional_c){

	double *RHS_i = get_col_Matrix_d(0, adjoint_data->RHS);

	int RHS_index; 
	double dI_by_dW_component; // The RHS components (I = Functional)

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
			dI_by_dW_component = cimag(functional_c(sim_c))/CX_STEP; 
			sol_coef_c->data[col_l] -= CX_STEP*I;

			// Load the value into the vector
			RHS_index = (int)s_vol_c->ind_dof+col_l;
			RHS_i[RHS_index] = dI_by_dW_component;
		}
	}
}



// static struct Multiarray_d* constructor_RHS_finite_diff(struct Optimization_Case *optimization_case){

// 	/*
// 	Compute the RHS ([dI/dW], where I is the objective function) using finite differences.
// 	Return the results in a multiarray

// 	Arguments:
// 		optimization_case = The data structure with the optimization material

// 	Return:
// 		Multiarray_d with the RHS components
// 	*/

// 	// Preprocessing:
// 	// Compute the size of the RHS vector

// 	struct Simulation *sim = optimization_case->sim; 

// 	int RHS_size_ex_0;  // Number of rows (extent 0) for RHS vector

// 	// Use the first volume as reference to get the P value
// 	// NOTE: Assumes global P is constant (need to take into account adaptation
// 	// eventually)
// 	struct Solver_Volume* s_vol = (struct Solver_Volume*)sim->volumes->first;

// 	int P = s_vol->p_ref,  // P value for all elements
// 		nv = (int)sim->n_v,  // number of volumes  
// 		neq = NEQ_EULER,  // number of equations
// 		d = DIM;  // Dimension

// 	RHS_size_ex_0 = nv * (int)(pow(P+1, d)) * neq;

	
// 	// =========================
// 	//       RHS = [dI/dW]
// 	// =========================

// 	// Compute the RHS vector using the complex step

// 	struct Multiarray_d *RHS = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){RHS_size_ex_0,1});
// 	double *RHS_i = get_col_Multiarray_d(0, RHS);

// 	int RHS_index; 
// 	double dI_by_dW_component;
// 	double I_0 = optimization_case->objective_function(sim);

// 	// Loop over the complex volumes, perturb the solution vector by a complex step and compute
// 	// each RHS component.
// 	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next){

// 		struct Solver_Volume *s_vol = (struct Solver_Volume*) curr;
// 		struct Multiarray_d* sol_coef = s_vol->sol_coef;

// 		// Perturb each solution vector component by the complex step and obtain the 
// 		// linearization terms.
// 		const ptrdiff_t n_col_l = compute_size(sol_coef->order,sol_coef->extents);
// 		for (int col_l = 0; col_l < n_col_l; col_l++) {
// 			sol_coef->data[col_l] += FINITE_DIFF_STEP;
// 			dI_by_dW_component = (optimization_case->objective_function(sim) - I_0)/FINITE_DIFF_STEP;
// 			sol_coef->data[col_l] -= FINITE_DIFF_STEP;

// 			// Load the value into the vector
// 			RHS_index = (int)s_vol->ind_dof+col_l;
// 			RHS_i[RHS_index] = dI_by_dW_component;
// 		}
// 	}

// 	return RHS;

// }


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


