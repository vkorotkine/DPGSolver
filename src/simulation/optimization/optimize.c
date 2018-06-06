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


#include "optimize.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

#include "face.h"
#include "volume.h"

#include "test_integration.h"

#include "definitions_adaptation.h"
#include "definitions_core.h"
#include "definitions_tol.h"

#include "macros.h"

#include "matrix.h"
#include "matrix_constructors.h"

#include "multiarray.h"
#include "complex_multiarray.h"
#include "multiarray_constructors.h"
#include "complex_multiarray_constructors.h"
#include "multiarray_math.h"

#include "intrusive.h"
#include "definitions_intrusive.h"

#include "geometry.h"
#include "geometry_parametric.h"
#include "solve_implicit.h"

#include "simulation.h"

#include "volume_solver.h"
#include "face_solver.h"
#include "test_complex_face_solver.h"
#include "test_complex_volume_solver.h"
#include "visualization.h"
#include "definitions_visualization.h"

#include "optimization_case.h"

#include "adjoint.h"
#include "sensitivities.h"
#include "optimization_minimizers.h"


// TODO: Read from the optimization.data file
#define CONST_L2_GRAD_EXIT 1E-8
#define MAX_NUM_DESIGN_ITERS 200

// Static function declarations ************************************************************************************* //


static struct Optimization_Case* setup_optimization(struct Simulation* sim);

static void copy_data_r_to_c_sim(struct Simulation* sim, struct Simulation* sim_c);

static struct Multiarray_d* compute_gradient(struct Optimization_Case* optimization_case);

static struct Multiarray_d* test_compute_brute_force_gradient(struct Optimization_Case* optimization_case);

// Interface functions ********************************************************************************************** //



void optimize(struct Simulation* sim){

	/*
	Optimize the the geometry to minimize a specified functional.

	Arguments:
		sim = The simulation data structure

	Return:
		- 
	*/

	// ================================
	//         Preprocessing
	// ================================
	struct Optimization_Case *optimization_case = setup_optimization(sim);
	copy_data_r_to_c_sim(optimization_case->sim, optimization_case->sim_c);


	// ================================
	//      Optimization Routine
	// ================================

	int design_iteration = 0;
	double L2_grad, objective_func_value;

	objective_func_value = optimization_case->objective_function(sim);
	printf(" INITIAL -> obj_func : %e  \n", objective_func_value);

	while(true){

		printf("------------------------------------\n");
		printf("Design Iteration : %d \n", design_iteration);fflush(stdout);

		// Setup and solve for the adjoint (Chi)
		setup_adjoint(optimization_case);
		solve_adjoint(optimization_case);


		// Compute dI_dXp and dR_dXp
		compute_sensitivities(optimization_case);
		printf("Computed Sensitivities \n");fflush(stdout);


		// Use the adjoint and sensitivities to find the gradient of the objective function
		optimization_case->grad_I = compute_gradient(optimization_case);  // free


		// Find the search direction, step length, and modify the geometry
		gradient_descent(optimization_case, design_iteration);


		// Solve the flow on the updated geometry
		solve_implicit(sim);


		// Compute optimization values to keep track and for the exit condition
		L2_grad = norm_Multiarray_d(optimization_case->grad_I, "L2");
		objective_func_value = optimization_case->objective_function(sim);
		printf(" L2_grad : %e   obj_func : %e  \n", L2_grad, objective_func_value);


		// Destruct allocated data structures:
		destruct_sensitivity_structures(optimization_case);
		destructor_Multiarray_d(optimization_case->grad_I);


		// Exit condition
		if (L2_grad < CONST_L2_GRAD_EXIT || design_iteration > MAX_NUM_DESIGN_ITERS)
			break;


		// Output the NURBS patch



		design_iteration++;

	}

	// ================================
	//         Postprocessing
	// ================================


	// Clear allocated structures:
	destructor_Optimization_Case(optimization_case);
	
	// Unused functions:
	UNUSED(test_compute_brute_force_gradient);

	/*
	double I_val = optimization_case->objective_function(optimization_case->sim);
	double complex I_val_c = optimization_case->objective_function_c(optimization_case->sim_c);
	
	printf("I   : %e \n", I_val);
	printf("I_c : %e + %e * i\n", creal(I_val_c), cimag(I_val_c));
	*/

}


// Static functions ************************************************************************************************* //


static void copy_data_r_to_c_sim(struct Simulation* sim, struct Simulation* sim_c){

	/*
	Copy the data from the real sim object to the complex sim object. This method will 
	transfer all the face solver and volume solver data into the complex sim object.

	NOTE: Since the mesh and same control file is read (no adaptation was done) we can 
		assume that the ordering of the volumes and faces in the sim and sim_c structure
		is the same.
	
	Arguments:
		sim = The sim data structure with real Volumes and Faces
		sim_c = The sim data structure with Volumes and Faces that hold complex data.

	Return:
		-
	*/

	// Copy the Volume_Solver data
	struct Intrusive_Link* curr   = sim->volumes->first;
	struct Intrusive_Link* curr_c = sim_c->volumes->first;

	while(true){

		struct Solver_Volume* s_vol 		= (struct Solver_Volume*) curr;
		struct Solver_Volume_c* s_vol_c 	= (struct Solver_Volume_c*) curr_c;
		
		copy_members_r_to_c_Solver_Volume((struct Solver_Volume_c*const)s_vol_c, 
			(const struct Solver_Volume*const)s_vol, sim);

		curr 	= curr->next;
		curr_c 	= curr_c->next;

		if(!curr || !curr_c)
			break;
	}

	// Copy the Face_Solver data
	curr   = sim->faces->first;
	curr_c = sim_c->faces->first;

	while(true){

		struct Solver_Face* s_face 		= (struct Solver_Face*) curr;
		struct Solver_Face_c* s_face_c 	= (struct Solver_Face_c*) curr_c;
		
		// TODO:
		// Needed in order to do the r_to_c conversion. Find a way to fix this issue
		s_face->nf_fc = (const struct const_Multiarray_d*)constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1,1});
		s_face_c->nf_fc = (const struct const_Multiarray_c*)constructor_empty_Multiarray_c('C',2,(ptrdiff_t[]){1,1});

		copy_members_r_to_c_Solver_Face((struct Solver_Face_c*const)s_face_c, 
			(const struct Solver_Face*const)s_face, sim);

		curr 	= curr->next;
		curr_c 	= curr_c->next;

		if(!curr || !curr_c)
			break;
	}

}


static struct Optimization_Case* setup_optimization(struct Simulation* sim){

	/*
	Setup the optimization case. In this method, first the constructor for the 
	Optimization_Case data structure will be called. Following this, a complex 
	sim object will be created. This object will be used for obtaining all linearizations
	using the complex step. The data from the real sim object will simply be copied into
	the complex counterpart between each design cycle.
	
	Arguments:
		sim = The simulation data structure
	
	Return:
		The Optimization_Case data structure holding the data needed for the optimization.

	*/

	struct Optimization_Case *optimization_case = constructor_Optimization_Case(sim);


	// Create the complex simulation object

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(sim->ctrl_name);
	
	const int* p_ref  = int_test_info->p_ref,
	         * ml_ref = int_test_info->ml;

	struct Simulation* sim_c = NULL;
	const char type_rc = 'c';

	bool ignore_static = false;
	int ml_max = ml_ref[1];

	UNUSED(ml_max);

	int ml_prev = ml_ref[0]-1,
    	p_prev  = p_ref[0]-1;

	int ml = ml_ref[0];
	int p = p_ref[0];

	const int adapt_type = int_test_info->adapt_type;

	structor_simulation(&sim_c, 'c', adapt_type, p, ml, p_prev, ml_prev, sim->ctrl_name, 
		type_rc, ignore_static);
	

	// Store pointers to the sim objects
	optimization_case->sim = sim;
	optimization_case->sim_c = sim_c;

	return optimization_case;

}


static struct Multiarray_d* compute_gradient(struct Optimization_Case* optimization_case){

	/*
	Compute the gradient of the objective function with respect to the design 
	variables. This is equal to 

		grad_I = dI_dXp - Chi_T * dR_dXp

	where Chi_T is the transpose of the adjoint.


	Arguments:
		optimization_case = The data structure holding the optimization information

	Return:
		A multiarray with the gradient of the objective function with respect to the 
		design variables. The multiarray is of dimension [1 x num_design_dofs] and 
		is in column major form.
	*/

	int num_design_dofs = (int)optimization_case->dI_dXp->extents[1];

	struct Multiarray_d *Chi_T, *grad_I, *Chi_T_dot_dR_dXp;
	grad_I = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, num_design_dofs});  // returned
	
	Chi_T = constructor_copy_Multiarray_d(optimization_case->Chi);  // destructed
	transpose_Multiarray_d(Chi_T, false);
	struct Matrix_d* Mat_Chi_T = constructor_move_Matrix_d_d('C',
				Chi_T->extents[0], Chi_T->extents[1], false, Chi_T->data); // destructed

	Chi_T_dot_dR_dXp = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, num_design_dofs});  // destructed

	mm_NN1C_Multiarray_d((const struct const_Matrix_d*const)Mat_Chi_T,
		(const struct const_Multiarray_d*const)optimization_case->dR_dXp ,Chi_T_dot_dR_dXp);

	for (int i = 0; i < num_design_dofs; i++){
		grad_I->data[i] = optimization_case->dI_dXp->data[i] - Chi_T_dot_dR_dXp->data[i];
	}

	// Destruct allocated structures
	destructor_Matrix_d(Mat_Chi_T);
	destructor_Multiarray_d(Chi_T_dot_dR_dXp);
	destructor_Multiarray_d(Chi_T);

	return grad_I;

}

static struct Multiarray_d* test_compute_brute_force_gradient(struct Optimization_Case* optimization_case){

	/*
	Compute the gradient using a brute force approach (perturb each design point and 
	run the flow solver). Use this gradient for testing purposes.

	Arguments:
		optimization_case = The data structure holding the optimization information

	Return:
		Multiarray of dimnesion [1 x num_design_dofs] corresponding to the gradient of
		the objective function with respect to the design variables.
	*/

	struct Simulation *sim = optimization_case->sim;

	int num_design_dofs = (int)optimization_case->dI_dXp->extents[1];
	struct Multiarray_d *grad_I = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, num_design_dofs});  // returned
	
	struct Multiarray_i* ctrl_pts_opt = optimization_case->geo_data.control_points_optimization;
	int *ctrl_pt_indeces = get_col_Multiarray_i(0, ctrl_pts_opt);
	int n_pts = (int)ctrl_pts_opt->extents[0];

	struct Multiarray_d* ctrl_pts_and_weights = optimization_case->geo_data.control_points_and_weights;
	int control_pt_index;

	double I_0 = optimization_case->objective_function(sim);
	int col_index = 0;

	for (int i = 0; i < n_pts; i++){
		// Loop over the design points

		for (int j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point
			// j = 1 (x degree of freeedom) and j = 2 (y degree of freedom)

			if (!get_col_Multiarray_i(j, ctrl_pts_opt)[i])
				continue;

			control_pt_index = ctrl_pt_indeces[i];
			get_col_Multiarray_d(j-1, ctrl_pts_and_weights)[control_pt_index] += FINITE_DIFF_STEP;


			update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)ctrl_pts_and_weights);
			set_up_solver_geometry(sim);
			solve_implicit(sim);
			grad_I->data[col_index] = (optimization_case->objective_function(sim) - I_0)/FINITE_DIFF_STEP;

			get_col_Multiarray_d(j-1, ctrl_pts_and_weights)[control_pt_index] -= FINITE_DIFF_STEP;

			col_index++;
			printf("col_index_brute_force : %d \n", col_index); fflush(stdout);
		}
	}

	// Update the geometry one last time with the original control point configurations
	update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)ctrl_pts_and_weights);
	set_up_solver_geometry(sim);
	solve_implicit(sim);

	return grad_I;

}





