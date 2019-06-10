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

#include "gradient.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "macros.h"
#include "definitions_tol.h"
#include "simulation.h"

#include "geometry.h"
#include "geometry_parametric.h"

#include "solve.h"
#include "solve_implicit.h"

#include "multiarray.h"
#include "multiarray_math.h"
#include "multiarray_constructors.h"
#include "matrix.h"
#include "matrix_constructors.h"

#include "optimization_case.h"
#include "functionals.h"


// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

struct Gradient_Data* constructor_Gradient_Data (struct Adjoint_Data *adjoint_data, 
	struct Sensitivity_Data *sensitivity_data){

	UNUSED(adjoint_data);

	struct Gradient_Data* gradient_data = calloc(1,sizeof *gradient_data); // returned

	int n_dof = (int) sensitivity_data->dI_dXp->ext_1;
	gradient_data->Gradient = constructor_empty_Matrix_d('C', 1, n_dof);

	return gradient_data;
}


void destructor_Gradient_Data(struct Gradient_Data* gradient_data){

	destructor_Matrix_d(gradient_data->Gradient);

	free((void*)gradient_data);

}


void compute_gradient(struct Gradient_Data *gradient_data, struct Adjoint_Data *adjoint_data, 
	struct Sensitivity_Data *sensitivity_data){

	int num_design_dofs = (int) gradient_data->Gradient->ext_1;

	// Setup the matrices for the gradient computation
	struct Matrix_d *Chi_T,
					*Chi_T_dot_dR_dXp;

	Chi_T = constructor_copy_Matrix_d(adjoint_data->Chi);  // free
	transpose_Matrix_d(Chi_T, false);
	
	Chi_T_dot_dR_dXp = constructor_empty_Matrix_d('C', 1, num_design_dofs);  // free

	// Compute {Chi}^T * [dR/dXp]
	mm_d('N', 'N', 1.0, 0.0, 
		(const struct const_Matrix_d*const)Chi_T,
		(const struct const_Matrix_d*const)sensitivity_data->dR_dXp,
		Chi_T_dot_dR_dXp);

	// Compute Gradient = [dI/dXp]^T - {Chi}^T * [dR/dXp]
	for (int i = 0; i < num_design_dofs; i++){
		gradient_data->Gradient->data[i] = sensitivity_data->dI_dXp->data[i] - Chi_T_dot_dR_dXp->data[i];
	}

	// Destruct allocated structures
	destructor_Matrix_d(Chi_T);
	destructor_Matrix_d(Chi_T_dot_dR_dXp);
}


void test_brute_force_gradient(struct Optimization_Case *optimization_case,
	functional_fptr functional, double *input_gradient){

	struct Simulation *sim = optimization_case->sim;

	int num_design_dofs = optimization_case->num_design_pts_dofs;
	struct Multiarray_d *grad_f = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, num_design_dofs});  // returned
	
	struct Multiarray_i* ctrl_pts_opt = optimization_case->geo_data.control_points_optimization;
	int *ctrl_pt_indeces = get_col_Multiarray_i(0, ctrl_pts_opt);
	int n_pts = (int)ctrl_pts_opt->extents[0];

	struct Multiarray_d* control_points = optimization_case->geo_data.control_points;
	int control_pt_index;

	double f_plus_h, f_min_h;  // Values of the functional for the central difference

	int col_index = 0;

	for (int i = 0; i < n_pts; i++){
		// Loop over the design points

		for (int j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point
			// j = 1 (x degree of freeedom) and j = 2 (y degree of freedom)

			if (!get_col_Multiarray_i(j, ctrl_pts_opt)[i])
				continue;

			control_pt_index = ctrl_pt_indeces[i];

			// f_plus_h
			get_col_Multiarray_d(j-1, control_points)[control_pt_index] += FINITE_DIFF_STEP;
			update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)control_points);
			set_up_solver_geometry(sim);
			solve_implicit(sim);
			f_plus_h = functional(sim);
			get_col_Multiarray_d(j-1, control_points)[control_pt_index] -= FINITE_DIFF_STEP;

			// f_min_h
			get_col_Multiarray_d(j-1, control_points)[control_pt_index] -= FINITE_DIFF_STEP;
			update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)control_points);
			set_up_solver_geometry(sim);
			solve_implicit(sim);
			f_min_h = functional(sim);
			get_col_Multiarray_d(j-1, control_points)[control_pt_index] += FINITE_DIFF_STEP;


			grad_f->data[col_index] = (f_plus_h - f_min_h)/(2.0*FINITE_DIFF_STEP);


			col_index++;
		}
	}

	// Update the geometry one last time with the original control point configurations
	//update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)control_points);
	//set_up_solver_geometry(sim);
	//solve_implicit(sim);

	printf("\n Gradient_Brute_Force: \n");
	for(int i = 0; i < num_design_dofs; i++)
		printf("%e ", grad_f->data[i]);

	printf("\n Gradient_Input: \n");
	for(int i = 0; i < num_design_dofs; i++)
		printf("%e ", input_gradient[i]);

	// Compute the L2 error
	int diff_vector_len = (int)compute_size(grad_f->order, grad_f->extents);

	struct Multiarray_d *difference_multiarray = constructor_empty_Multiarray_d(
		'C',2,(ptrdiff_t[]){1, diff_vector_len});  // free

	for (int i = 0; i < diff_vector_len; i++)
		difference_multiarray->data[i] = grad_f->data[i] - input_gradient[i];

	printf("\n");
	printf("Difference L2 Norm: %e \n",  norm_Multiarray_d(difference_multiarray, "L2"));


}

