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

#include "optimization_minimizers.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "macros.h"

#include "optimization_case.h"
#include "simulation.h"

#include "multiarray.h"

#include "matrix.h"
#include "matrix_constructors.h"
#include "matrix_math.h"

#include "geometry.h"
#include "geometry_parametric.h"

#include "solve_implicit.h"


// TODO: Read this from the optimization.data file
#define CONST_GRADIENT_DESCENT_STEP_SIZE_INIT 1E-3
#define CONST_BFGS_STEP_SIZE_INIT 1.0
#define CONST_WOLFE_CONDITION_C 1E-4
#define CONST_WOLFE_CONDITION_RHO 1E-1
#define CONST_WOLFE_CONDITION_ALPHA_MIN 1E-9
#define MAX_NORM_P 1E-2
#define NUM_SMOOTHING_ITER_BFGS 4

// Static function declarations ************************************************************************************* //

static void output_NURBS_patch_information(struct Optimization_Case* optimization_case);

static double backtracking_line_search(struct Optimization_Case* optimization_case, 
	double* p_k, double* grad_f_k, double alpha);

static void update_design_points(struct Optimization_Case* optimization_case, double *p_k, double alpha);

// Interface functions ********************************************************************************************** //


void preprocessor_minimizer(struct Optimization_Case* optimization_case){

	/*
	This function performs any preprocessing needed for the optimization minimizers.
	The memory for data structures needed by the optimizers is set up in this method. 
	Note that the data structures for all the different 

	Arguments:
		optimization_case = The data structure holding the optimization information
	*/

	int num_design_pts_dofs = optimization_case->num_design_pts_dofs;

	// Smoothing Matrix:
	// Create a laplacian smoothing matrix with epsilon = 0.6
	optimization_case->Laplacian_Matrix = constructor_empty_Matrix_d('C', num_design_pts_dofs, num_design_pts_dofs);
	
	set_to_value_Matrix_d(optimization_case->Laplacian_Matrix, 0.0);
	double eps_val = 0.6;

	for (int i = 0; i < num_design_pts_dofs; i++){

		get_col_Matrix_d(i, optimization_case->Laplacian_Matrix)[i] = -1.0 * eps_val;
		if (i == 0){
			// First row
			get_col_Matrix_d(i+1, optimization_case->Laplacian_Matrix)[i] = 0.5 * eps_val;
		} else if (i == num_design_pts_dofs-1){
			// Last row
			get_col_Matrix_d(i-1, optimization_case->Laplacian_Matrix)[i] = 0.5 * eps_val;
		} else{
			// Middle row
			get_col_Matrix_d(i+1, optimization_case->Laplacian_Matrix)[i] = 0.5 * eps_val;
			get_col_Matrix_d(i-1, optimization_case->Laplacian_Matrix)[i] = 0.5 * eps_val;
		}
	}

	// BFGS Data Structures

	optimization_case->B_k_inv = constructor_empty_Matrix_d('C', num_design_pts_dofs, num_design_pts_dofs);
	optimization_case->s_k = constructor_empty_Matrix_d('C', num_design_pts_dofs, 1);
	optimization_case->grad_f_k = constructor_empty_Matrix_d('C', num_design_pts_dofs, 1);

	// - Set the B_kMin1_inv matrix to be an identity matrix first
	set_to_value_Matrix_d(optimization_case->B_k_inv, 0.0);
	for(int i = 0; i < num_design_pts_dofs; i++)
		get_col_Matrix_d(i, optimization_case->B_k_inv)[i] = 1.0;
}


void postprocessor_minimizer(struct Optimization_Case* optimization_case){

	/*
	This function performs any postprocessing needed for the optimization minimizers.
	Data allocated for use by the minimizers in the preprocessor is freed here.

	Arguments:
		optimization_case = The data structure holding the optimization information
	*/

	// Smoothing Matrix

	destructor_Matrix_d(optimization_case->Laplacian_Matrix);

	// BFGS Data Structures

	destructor_Matrix_d(optimization_case->B_k_inv);
	destructor_Matrix_d(optimization_case->s_k);
	destructor_Matrix_d(optimization_case->grad_f_k);

}


void BFGS_minimizer(struct Optimization_Case *optimization_case, int design_iteration){

	/*
	Use the BFGS method to compute the search direction for the optimization. Then, use
	the search direction to compute the step length using the backtracking 
	algorithm.

	Arguments:
		optimization_case = The optimization case data structure. Holds the information
			about the gradient of the objective.
		design_iteration = The current design iteration.

	Return:
		-
	*/

	int num_design_pts_dofs = optimization_case->num_design_pts_dofs;

	struct Matrix_d *p_k, *p_k_smooth;  // search direction (and smoothing search direction)
	p_k 		= constructor_empty_Matrix_d('C',num_design_pts_dofs,1);  // free
	p_k_smooth 	= constructor_empty_Matrix_d('C',num_design_pts_dofs,1);  // free

	if(design_iteration > 0){
		
		// Compute the Hessian approximation using the information for the previous
		// design step (k) to find the information for this design step (k+1)

		struct Matrix_d *B_k_inv 	= optimization_case->B_k_inv;
		struct Matrix_d *s_k 		= optimization_case->s_k;
		struct Matrix_d *grad_f_k 	= optimization_case->grad_f_k;

		double skT_dot_yk, ykT_dot_BKInv_dot_yk; 

		struct Matrix_d *y_k 					= constructor_empty_Matrix_d('C',num_design_pts_dofs,1);  // free
		struct Matrix_d *BKInv_dot_yk 			= constructor_empty_Matrix_d('C',num_design_pts_dofs,1);  // free
		struct Matrix_d *yk_dot_skT 			= constructor_empty_Matrix_d('C',num_design_pts_dofs, num_design_pts_dofs);  // free
		struct Matrix_d *BKInv_dot_yk_dot_skT 	= constructor_empty_Matrix_d('C',num_design_pts_dofs, num_design_pts_dofs);  // free
		struct Matrix_d *sk_dot_ykT 			= constructor_empty_Matrix_d('C',num_design_pts_dofs, num_design_pts_dofs);  // free
		struct Matrix_d *sk_dot_ykT_dot_BKInv 	= constructor_empty_Matrix_d('C',num_design_pts_dofs, num_design_pts_dofs);  // free
		struct Matrix_d *sk_dot_skT 			= constructor_empty_Matrix_d('C',num_design_pts_dofs, num_design_pts_dofs);  // free
		
		// Compute y_k
		// y_k = grad_f_kPlus1 - grad_f_k
		//  - NOTE: optimization_case->grad_I holds the gradient for this design step (k+1)
		// 		so it is equivalent to grad_f_kPlus1
		for (int i = 0; i < num_design_pts_dofs; i++)
			y_k->data[i] = optimization_case->grad_I->data[i] - grad_f_k->data[i];


		// Compute skT_dot_yk
		skT_dot_yk = 0;	
		for(int i = 0; i < num_design_pts_dofs; i++){
			skT_dot_yk += s_k->data[i] * y_k->data[i];
		}


		// Compute ykT_dot_BKInv_dot_yk
		mm_d('N', 'N', 1.0, 0.0, 
			(const struct const_Matrix_d*const)B_k_inv, 
			(const struct const_Matrix_d*const)y_k, 
			BKInv_dot_yk);

		ykT_dot_BKInv_dot_yk = 0;
		for(int i = 0; i < num_design_pts_dofs; i++){
			ykT_dot_BKInv_dot_yk += y_k->data[i] * BKInv_dot_yk->data[i];
		}


		// Compute sk_dot_skT
		mm_d('N', 'T', 1.0, 0.0, 
			(const struct const_Matrix_d*const)s_k, 
			(const struct const_Matrix_d*const)s_k, 
			sk_dot_skT);


		// Compute BKInv_dot_yk_dot_skT
		mm_d('N', 'T', 1.0, 0.0, 
			(const struct const_Matrix_d*const)y_k, 
			(const struct const_Matrix_d*const)s_k, 
			yk_dot_skT);

		mm_d('N', 'N', 1.0, 0.0, 
			(const struct const_Matrix_d*const)B_k_inv, 
			(const struct const_Matrix_d*const)yk_dot_skT, 
			BKInv_dot_yk_dot_skT);


		// Compute sk_dot_ykT_dot_BKInv
		mm_d('N', 'T', 1.0, 0.0, 
			(const struct const_Matrix_d*const)s_k, 
			(const struct const_Matrix_d*const)y_k, 
			sk_dot_ykT);

		mm_d('N', 'N', 1.0, 0.0, 
			(const struct const_Matrix_d*const)sk_dot_ykT, 
			(const struct const_Matrix_d*const)B_k_inv, 
			sk_dot_ykT_dot_BKInv);


		// Compute B_inv_kPlus1. Place these values into B_k_inv since these new
		// values will form the inverse hessian for this design step coming up
		for (int i = 0; i < num_design_pts_dofs * num_design_pts_dofs; i++){

			B_k_inv->data[i] = B_k_inv->data[i]
				+ ((skT_dot_yk + ykT_dot_BKInv_dot_yk)/(skT_dot_yk*skT_dot_yk)) * sk_dot_skT->data[i] 
				- (1./skT_dot_yk) * (BKInv_dot_yk_dot_skT->data[i] + sk_dot_ykT_dot_BKInv->data[i]);
		}

		// Free allocated structures
		destructor_Matrix_d(y_k);
		destructor_Matrix_d(BKInv_dot_yk);
		destructor_Matrix_d(yk_dot_skT);
		destructor_Matrix_d(BKInv_dot_yk_dot_skT);
		destructor_Matrix_d(sk_dot_ykT);
		destructor_Matrix_d(sk_dot_ykT_dot_BKInv);
		destructor_Matrix_d(sk_dot_skT);

	}

	// Now that the Hessian for this design step has been found, store the data 
	// for this kth step so that we can use it to get the next steps approximation

	// Save grad_f_k so that B_kPlus1_inv can be found (next design iteration)
	for (int i = 0; i < num_design_pts_dofs; i++)
		optimization_case->grad_f_k->data[i] = optimization_case->grad_I->data[i];


	// p_k = - B_k_inv * grad_f_k;
	mm_d('N', 'N', -1.0, 0.0,
		(const struct const_Matrix_d*const)optimization_case->B_k_inv,
		(const struct const_Matrix_d*const)optimization_case->grad_f_k,
		p_k);

	// TEMPORARY ADDITIONs:

	// 1) Apply the smoothing
	for (int i = 0; i < NUM_SMOOTHING_ITER_BFGS; i++){
		mm_d('N', 'N', 1.0, 0.0,
			(const struct const_Matrix_d*const)optimization_case->Laplacian_Matrix,
			(const struct const_Matrix_d*const)p_k,
			p_k_smooth);

		for (int j = 0; j < num_design_pts_dofs; j++)
			p_k->data[j] = p_k_smooth->data[j];

	}

	// 2) Make the max norm (euclidean) of the P vector MAX_NORM_P
	double norm_P = 0.0;
	for (int i = 0; i < num_design_pts_dofs; i++){
		norm_P += p_k->data[i] * p_k->data[i];
	}
	norm_P = sqrt(norm_P);
	
	if (norm_P > MAX_NORM_P){
		printf("Scale -> norm_P : %e \n", norm_P);
		for (int i = 0; i < num_design_pts_dofs; i++){
			p_k->data[i] = (MAX_NORM_P/norm_P) * p_k->data[i];
		}	
	}


	// Perform the back tracking line search
	double alpha_k = backtracking_line_search(optimization_case, p_k->data, optimization_case->grad_f_k->data, 
		CONST_BFGS_STEP_SIZE_INIT);


	// Save s_k = alpha_k * p_k so that B_kPlus1_inv can be found (next design iteration)
	for (int i = 0; i < num_design_pts_dofs; i++)
		optimization_case->s_k->data[i] = alpha_k * p_k->data[i];


	// Destruct allocated data structures:
	destructor_Matrix_d(p_k);
	destructor_Matrix_d(p_k_smooth);

}



void gradient_descent(struct Optimization_Case *optimization_case, int design_iteration){

	/*
	Traverse in the negative gradient direction to minimize the functional. This function
	will move a small step length in the negative gradient direction. 

	Arguments:
		optimization_case = The optimization case data structure. Holds the information
			about the gradient of the objective.
		design_iteration = The current design iteration.

	Return:
		-
	*/

	UNUSED(design_iteration);

	// Compute the search direction. For gradient descent, this is equal to the negative of grad_I

	int num_design_pts_dofs = optimization_case->num_design_pts_dofs;
	struct Multiarray_d *p_k;  // search direction
	p_k = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_design_pts_dofs,1});  // free

	for (int i = 0; i < num_design_pts_dofs; i++)
		p_k->data[i] = -1.0 * optimization_case->grad_I->data[i];

	// Perform the back tracking line search
	backtracking_line_search(optimization_case, p_k->data, optimization_case->grad_I->data, 
		CONST_GRADIENT_DESCENT_STEP_SIZE_INIT);

	// Destruct allocated data structures:
	destructor_Multiarray_d(p_k);

	printf(" Completed Gradient Descent \n");  // Monitor Progress

}


// Static functions ************************************************************************************************* //


static void update_design_points(struct Optimization_Case* optimization_case, double *p_k, double alpha){

	/*
	Traverse along the search direction p_k by the step length amount alpha. Modify the 
	control points.
		
		x_{k+1} = x_k + alpha * p_k

	x_{k+1} = Modified design control points
	x_k = current design control points


	Arguments:
		optimization_case = The data structure holding the optimization information
		p_k = The search direction. Should have num_design_pts_dofs values.
		alpha = The step length

	Return:
		-
	*/

	struct Multiarray_d *ctrl_pts_and_weights = optimization_case->geo_data.control_points_and_weights;
	struct Multiarray_i *ctrl_pts_optimization = optimization_case->geo_data.control_points_optimization;
	int n_pts_optimization = (int)ctrl_pts_optimization->extents[0];

	int control_pt_index;
	int p_index = 0;

	for (int i = 0; i < n_pts_optimization; i++){
		// Loop over the design optimization points

		for (int j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point
			// j = 1 (x degree of freeedom) and j = 2 (y degree of freedom)

			if (!get_col_Multiarray_i(j, ctrl_pts_optimization)[i])
				continue;

			// Modify the control point position
			control_pt_index = get_col_Multiarray_i(0, ctrl_pts_optimization)[i];
			get_col_Multiarray_d(j-1, ctrl_pts_and_weights)[control_pt_index] += alpha * p_k[p_index++];
		}
	}

	// Monitor the progress
	output_NURBS_patch_information(optimization_case);

}


static double backtracking_line_search(struct Optimization_Case* optimization_case, double* p_k, double* grad_f_k,
	double alpha){

	/*
	Perform the back tracking line search and ensure that the first Wolfe condition
	is satisfied. This function will exit once the optimal alpha (step length) is computed.
	Thus, upon exiting the function, the geometry will have been updated to the desired
	shape for the given design iteration.

	Arguments:
		optimization_case = The data structure holding the optimization information.
		p_k	= The search direction for this iteration (kth iteration) as an array of 
			doubles with num_design_pts_dofs values.
		grad_f_k = The gradient of the objective function for this iteration (kth iteration).
		alpha = The initial step length value

	Return:
		Alpha value used as the step length
	*/

	struct Simulation *sim = optimization_case->sim;

	int num_design_pts_dofs = optimization_case->num_design_pts_dofs;

	// Compute the dot product of the search direction and gradient
	double p_k_dot_grad_f_k = 0; 
	for (int i = 0; i < num_design_pts_dofs; i++)
		p_k_dot_grad_f_k += p_k[i]*grad_f_k[i];


	// The value of the objective function (f) at x_k and x_k + alpha*p_k
	double f_x_k, f_x_k_plus_alpha_p_k;
	f_x_k = optimization_case->objective_function(sim);

	while (true){

		// Move along the search direction. Update the geometry and recompute the metric terms
		update_design_points(optimization_case, p_k, alpha);
		update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)optimization_case->geo_data.control_points_and_weights);
		set_up_solver_geometry(sim);

		// Solve the flow on the updated geometry
		solve_implicit(sim);

		// Check Wolfe Condition
		f_x_k_plus_alpha_p_k = optimization_case->objective_function(sim);
		if (f_x_k_plus_alpha_p_k <= (f_x_k + CONST_WOLFE_CONDITION_C*alpha*p_k_dot_grad_f_k))
			break;

		// If the step length was below a constant minimum value stop the line search
		if (alpha < CONST_WOLFE_CONDITION_ALPHA_MIN){
			printf("\n\n EXITING BACKTRACK -> CONST_WOLFE_CONDITION_ALPHA_MIN \n\n");
			break;
		}

		// Step length not valid (and not too small), so revert the geometry design points to their previous location
		update_design_points(optimization_case, p_k, -1.0*alpha);

		// Decrease step length by factor rho
		printf("\n\n (%e, %e) Alpha : %e -> Alpha : %e \n\n", f_x_k_plus_alpha_p_k, 
			f_x_k + CONST_WOLFE_CONDITION_C*alpha*p_k_dot_grad_f_k, alpha, alpha*CONST_WOLFE_CONDITION_RHO);  // monitor progress
		alpha *= CONST_WOLFE_CONDITION_RHO;

	}

	printf("\n\n Completed backtracking_line_search\n\n");  // monitor progress

	return alpha;

}

static void output_NURBS_patch_information(struct Optimization_Case* optimization_case){

	struct Simulation *sim = optimization_case->sim;

	char f_name[4*STRLEN_MAX] = { 0, };
	sprintf(f_name,"%s%c%s%c%s", sim->pde_name,'/',sim->pde_spec,'/',"NURBS_Patch");

	char output_name[STRLEN_MAX] = { 0, };
	strcpy(output_name,"../output/");
	strcat(output_name,"paraview/");
	strcat(output_name,f_name);

	FILE* fp;
	if ((fp = fopen(output_name,"w")) == NULL)
		printf("Error: File %s did not open.\n", output_name), exit(1);

	// Print the patch information
	fprintf(fp, "P(xi_order) %d\n", optimization_case->geo_data.P);
	fprintf(fp, "Q(eta_order) %d\n", optimization_case->geo_data.Q);
	fprintf(fp, "\n");

	fprintf(fp, "knots_xi %d\n", (int)optimization_case->geo_data.knots_xi->extents[0]);
	for (int i = 0; i < (int)optimization_case->geo_data.knots_xi->extents[0]; i++){
		fprintf(fp, "%.14e\n", optimization_case->geo_data.knots_xi->data[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "knots_eta %d\n", (int)optimization_case->geo_data.knots_eta->extents[0]);
	for (int i = 0; i < (int)optimization_case->geo_data.knots_eta->extents[0]; i++){
		fprintf(fp, "%.14e\n", optimization_case->geo_data.knots_eta->data[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Control_Point_Data %d \n", (int)optimization_case->geo_data.control_points_and_weights->extents[0]);
	for (int i = 0; i < (int)optimization_case->geo_data.control_points_and_weights->extents[0]; i++){
		fprintf(fp, "%.14e %.14e %.14e\n", 
			get_col_Multiarray_d(0, optimization_case->geo_data.control_points_and_weights)[i],
			get_col_Multiarray_d(1, optimization_case->geo_data.control_points_and_weights)[i],
			get_col_Multiarray_d(2, optimization_case->geo_data.control_points_and_weights)[i]
			);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Control_Point_Connectivity %d %d\n", 
		(int)optimization_case->geo_data.control_point_connectivity->extents[0],
		(int)optimization_case->geo_data.control_point_connectivity->extents[1]);
	for (int i = 0; i < (int)optimization_case->geo_data.control_point_connectivity->extents[0]; i++){
		for (int j = 0; j < (int)optimization_case->geo_data.control_point_connectivity->extents[1]; j++){
			fprintf(fp, "%d ", get_row_Multiarray_i(i, optimization_case->geo_data.control_point_connectivity)[j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);

}


