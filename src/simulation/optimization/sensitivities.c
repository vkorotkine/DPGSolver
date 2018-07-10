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


#include "sensitivities.h"

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
#include "intrusive.h"
#include "volume_solver.h"

#include "geometry.h"
#include "geometry_parametric.h"

#include "multiarray.h"
#include "multiarray_constructors.h"
#include "matrix.h"
#include "matrix_constructors.h"

#include "solve.h"
#include "solve_implicit.h"
#include "compute_grad_coef_dg.h"
#include "compute_volume_rlhs_dg.h"
#include "compute_face_rlhs_dg.h"
#include "computational_elements.h"
#include "definitions_intrusive.h"
#include "test_case.h"

#include "optimization_case.h"
#include "functionals.h"


// Static function declarations ************************************************************************************* //


/** \brief 	Compute the sensitivities of the functional (I) and residual (R) with respect to the 
 *	design variables (Xp). This method computes the components of the dI_dXp structure 
 *	and dR_dXp structures using finite differences. 
 *
 *	Ordering of sensitivities
 *		- If there are n control points on the NURBS patch to use for the optimization, 
 *			each point may be moved in either the x or y direction (only 2D patches 
 *			considered for now). Information for the freedom of each point is stored in 
 *			control_points_optimization in optimization_case->geo_data. 
 *		- The order is found as follows. For each design point, the partial with respect to the
 *			x degree of freedom and then the y degree of freedom are listed. This is continued
 *			for each control point (the ordering of the control points is given by the row ordering
 * 			of the points in control_points_optimization)
 *		- example:
 *		 	if control_points_optimization = [[1, 0, 1], [2, 1, 1]] then:
 *
 *			 	The first col of dI_dXp is the partial with respect to Xp1's y degree of freedom, 
 *			 	the second col is the partial with respect to Xp2's x degree of freedom and 
 *			 	the third col is the partial with respect to Xp2's y degree of freedom.
 */
static void set_dI_dXp_and_dR_dXp_finite_diff(
	struct Sensitivity_Data *sensitivity_data, ///< The sensitivity_data data structure to load data into
	struct Optimization_Case *optimization_case, ///< Standard. The optimization_case data structure
	functional_fptr functional ///< The functional function pointer (real version)
	);


/** \brief 	Compute the sensitivities of the functional (I) and residual (R) with respect to the 
 *	design variables (Xp). This method computes the components of the dI_dXp structure 
 *	and dR_dXp structures using the complex step. 
 *
 *	Ordering of sensitivities
 *		- If there are n control points on the NURBS patch to use for the optimization, 
 *			each point may be moved in either the x or y direction (only 2D patches 
 *			considered for now). Information for the freedom of each point is stored in 
 *			control_points_optimization in optimization_case->geo_data. 
 *		- The order is found as follows. For each design point, the partial with respect to the
 *			x degree of freedom and then the y degree of freedom are listed. This is continued
 *			for each control point (the ordering of the control points is given by the row ordering
 * 			of the points in control_points_optimization)
 *		- example:
 *		 	if control_points_optimization = [[1, 0, 1], [2, 1, 1]] then:
 *
 *			 	The first col of dI_dXp is the partial with respect to Xp1's y degree of freedom, 
 *			 	the second col is the partial with respect to Xp2's x degree of freedom and 
 *			 	the third col is the partial with respect to Xp2's y degree of freedom.
 */
static void set_dI_dXp_and_dR_dXp_cmplx_step(
	struct Sensitivity_Data *sensitivity_data, ///< The sensitivity_data data structure to load data into
	struct Optimization_Case *optimization_case, ///< Standard. The optimization_case data structure
	functional_fptr_c functional_c ///< The functional function pointer (real version)
	);


/** \brief 	Compute the residual vector. Do this by calling compute_rlhs and copying the 
 *	negative of the RHS vector in the implicit solve.
 *
 * \return 	A multiarray of dimension [num_res_eq x 1] with the value of the resiudal 
 * 	vector.
 */
static struct Multiarray_d* constructor_residual(
	struct Sensitivity_Data *sensitivity_data, ///< Standard. The sensitivity_data data structure
	struct Optimization_Case *optimization_case ///< Standard. The optimization_case data structure
	);


static struct Multiarray_c* constructor_residual_c(
	struct Sensitivity_Data *sensitivity_data, ///< Standard. The sensitivity_data data structure
	struct Optimization_Case *optimization_case ///< Standard. The optimization_case data structure
	);



// Interface functions ********************************************************************************************** //


struct Sensitivity_Data* constructor_Sensitivity_Data (struct Optimization_Case *optimization_case){

	struct Sensitivity_Data* sensitivity_data = calloc(1,sizeof *sensitivity_data); // returned


	struct Simulation *sim = optimization_case->sim;
	
	int num_design_pt_dofs = optimization_case->num_design_pts_dofs; // The total number of design point dofs
	
	int num_res_eqs = 0; // The total number of residual equations (and flow variables)
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;
		num_res_eqs += (int)compute_size(s_vol->sol_coef->order,s_vol->sol_coef->extents);
	}

	// Allocate the memory for the dI_dXp and dR_dXp structures
	sensitivity_data->dR_dXp = constructor_empty_Matrix_d('C', num_res_eqs, num_design_pt_dofs);
	sensitivity_data->dI_dXp = constructor_empty_Matrix_d('C', 1, num_design_pt_dofs);

	return sensitivity_data;

}


void destructor_Sensitivity_Data(struct Sensitivity_Data* sensitivity_data){

	destructor_Matrix_d(sensitivity_data->dR_dXp);
	destructor_Matrix_d(sensitivity_data->dI_dXp);

	free((void*)sensitivity_data);

}


void compute_sensitivities(struct Sensitivity_Data* sensitivity_data, struct Optimization_Case *optimization_case, 
	functional_fptr functional, functional_fptr_c functional_c){


	UNUSED(functional);
	UNUSED(functional_c);

	UNUSED(set_dI_dXp_and_dR_dXp_finite_diff);
	UNUSED(set_dI_dXp_and_dR_dXp_cmplx_step);

	// set_dI_dXp_and_dR_dXp_finite_diff(sensitivity_data, optimization_case, functional);

	// struct Matrix_d *dR_dXp_finite_diff = constructor_copy_Matrix_d(sensitivity_data->dR_dXp); // free
	// struct Matrix_d *dI_dXp_finite_diff = constructor_copy_Matrix_d(sensitivity_data->dI_dXp); // free

	// //printf("dR_dXp_finite_diff: \n");
	//print_Matrix_d(dR_dXp_finite_diff);

	// printf("dI_dXp_finite_diff: \n");
	// print_Matrix_d(dI_dXp_finite_diff);

	// destructor_Matrix_d(dR_dXp_finite_diff);
	// destructor_Matrix_d(dI_dXp_finite_diff);

	set_dI_dXp_and_dR_dXp_cmplx_step(sensitivity_data, optimization_case, functional_c);
	//set_dI_dXp_and_dR_dXp_finite_diff(sensitivity_data, optimization_case, functional);

	// printf("dI_dXp\n");
	// print_Matrix_d(sensitivity_data->dI_dXp);

	//printf("dR_dXp complex step\n");
	//print_Matrix_d(sensitivity_data->dR_dXp);

	// exit(0);
}


// Static functions ************************************************************************************************* //


static void set_dI_dXp_and_dR_dXp_finite_diff(struct Sensitivity_Data *sensitivity_data, struct Optimization_Case *optimization_case,
	functional_fptr functional){

	struct Simulation *sim = optimization_case->sim;
	
	struct Multiarray_i* ctrl_pts_opt = optimization_case->geo_data.control_points_optimization;
	int *ctrl_pt_indeces = get_col_Multiarray_i(0, ctrl_pts_opt);
	int n_pts = (int)ctrl_pts_opt->extents[0];


 	struct Matrix_d *dI_dXp = sensitivity_data->dI_dXp;
	double *dI_dXp_i = dI_dXp->data;

	struct Matrix_d *dR_dXp = sensitivity_data->dR_dXp;
	int num_res_eqs = (int) dR_dXp->ext_0;

	// =================================================
	//     Compute the dI_dXp and dR_dXp components
	// =================================================


	// Compute initial I_0 (objective func) and residual (R_0) values
	double I_0 					= functional(sim);
	struct Multiarray_d *R_0 	= constructor_residual(sensitivity_data, optimization_case);  // free

	struct Multiarray_d* control_points = optimization_case->geo_data.control_points;
	int control_pt_index;

	int col_index = 0;

	for (int i = 0; i < n_pts; i++){
		// Loop over the design points

		for (int j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point
			// j = 1 (x degree of freeedom) and j = 2 (y degree of freedom)

			if (!get_col_Multiarray_i(j, ctrl_pts_opt)[i])
				continue;

			control_pt_index = ctrl_pt_indeces[i];
			get_col_Multiarray_d(j-1, control_points)[control_pt_index] += FINITE_DIFF_STEP;


			update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)control_points);
			set_up_solver_geometry(sim);
			dI_dXp_i[col_index] = (functional(sim) - I_0)/FINITE_DIFF_STEP;
			
			struct Multiarray_d *R = constructor_residual(sensitivity_data, optimization_case);  // free
			
			for (int k = 0; k < num_res_eqs; k++)
				get_col_Matrix_d(col_index, dR_dXp)[k] = 
					(1./FINITE_DIFF_STEP)*(get_col_Multiarray_d(0, R)[k] - get_col_Multiarray_d(0, R_0)[k]);
			
			destructor_Multiarray_d(R);
			

			get_col_Multiarray_d(j-1, control_points)[control_pt_index] -= FINITE_DIFF_STEP;

			col_index++;
		}
	}

	// Update the geometry one last time with the original control point configurations
	update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)control_points);
	set_up_solver_geometry(sim);

	// Free allocated structures
	destructor_Multiarray_d(R_0);
}


static struct Multiarray_d* constructor_residual(struct Sensitivity_Data* sensitivity_data, 
	struct Optimization_Case *optimization_case){

	// MSB: TODO: Set residual without calling PETSc constructor. Instead, compute the
	// residual for each equation and load it manually here into the mutliarray (so that 
	// the complex version can be easily done)

	struct Simulation* sim = optimization_case->sim;
	struct Solver_Storage_Implicit* ssi = constructor_Solver_Storage_Implicit(sim); // destructed
	compute_rlhs_optimization(sim, ssi);

	// Set Residual Multiarray
	int num_res_eqs = (int)sensitivity_data->dR_dXp->ext_0;
	struct Multiarray_d* Residual = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_res_eqs, 1});
	double *Residual_i = get_col_Multiarray_d(0, Residual);
	

	// TODO: Find a way to find residual without getting the LHS also (and filling it into PETSc)


	// Copy values into Residual Multiarray (recall RHS in implicit solve is the negative of the residual)
	int *ix;
	ix = malloc((unsigned int)num_res_eqs* sizeof *ix);  // free

	for (int i = 0; i < num_res_eqs; i++)
		ix[i] = i;

	VecGetValues(ssi->b, num_res_eqs, ix, Residual_i);

	for (int i = 0; i < num_res_eqs; i++)
		Residual_i[i] *= -1.;

	// Free allocated structures
	free(ix);
	destructor_Solver_Storage_Implicit(ssi);


	return Residual;

}



static struct Multiarray_c* constructor_residual_c(struct Sensitivity_Data* sensitivity_data, 
	struct Optimization_Case *optimization_case){

	// MSB: TODO: Set residual without calling PETSc constructor. Instead, compute the
	// residual for each equation and load it manually here into the mutliarray (so that 
	// the complex version can be easily done)

	struct Simulation* sim_c = optimization_case->sim_c;

	((struct Test_Case_c*)sim_c->test_case_rc->tc)->solver_method_curr = 'e';
	constructor_derived_Elements(sim_c,IL_ELEMENT_SOLVER_DG); // destructed
	constructor_derived_computational_elements_c(sim_c,IL_SOLVER_DG); // destructed
	
	initialize_zero_memory_volumes_c(sim_c->volumes);
	compute_grad_coef_dg_c(sim_c,sim_c->volumes,sim_c->faces);
	compute_volume_rlhs_dg_c(sim_c,NULL,sim_c->volumes);
	compute_face_rlhs_dg_c(sim_c,NULL,sim_c->faces);
	
	destructor_derived_Elements(sim_c,IL_ELEMENT_SOLVER);
	destructor_derived_computational_elements_c(sim_c,IL_SOLVER);
	((struct Test_Case_c*)sim_c->test_case_rc->tc)->solver_method_curr = 'i';

	// Set Residual Multiarray
	int num_res_eqs = (int)sensitivity_data->dR_dXp->ext_0;
	struct Multiarray_c* Residual = constructor_empty_Multiarray_c('C',2,(ptrdiff_t[]){num_res_eqs, 1});
	double complex *Residual_i = get_col_Multiarray_c(0, Residual);
	

	for (struct Intrusive_Link* curr_c = sim_c->volumes->first; curr_c; curr_c = curr_c->next) {
		struct Solver_Volume_c* s_vol_c = (struct Solver_Volume_c*) curr_c;
		
		struct Multiarray_c* sol_coef_c = s_vol_c->sol_coef;
		struct Multiarray_c* rhs_c = s_vol_c->rhs;

		ptrdiff_t ext_0 = compute_size(sol_coef_c->order,sol_coef_c->extents);

		for (int i = 0; i < ext_0; i++){
			Residual_i[(int)s_vol_c->ind_dof + i] = rhs_c->data[i];
		}
	}

	return Residual;
}




static void set_dI_dXp_and_dR_dXp_cmplx_step(struct Sensitivity_Data *sensitivity_data, struct Optimization_Case *optimization_case,
	functional_fptr_c functional_c){

	struct Simulation *sim_c = optimization_case->sim_c;
	
	struct Multiarray_i* ctrl_pts_opt = optimization_case->geo_data.control_points_optimization;
	int *ctrl_pt_indeces = get_col_Multiarray_i(0, ctrl_pts_opt);
	int n_pts = (int)ctrl_pts_opt->extents[0];


 	struct Matrix_d *dI_dXp = sensitivity_data->dI_dXp;
	double *dI_dXp_i = dI_dXp->data;
	UNUSED(dI_dXp_i);
	UNUSED(functional_c);

	struct Matrix_d *dR_dXp = sensitivity_data->dR_dXp;
	int num_res_eqs = (int) dR_dXp->ext_0;

	// =================================================
	//     Compute the dI_dXp and dR_dXp components
	// =================================================

	struct Multiarray_c* control_points_c = optimization_case->geo_data.control_points_c;
	int control_pt_index;

	int col_index = 0;

	for (int i = 0; i < n_pts; i++){
		// Loop over the design points

		for (int j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point
			// j = 1 (x degree of freeedom) and j = 2 (y degree of freedom)

			if (!get_col_Multiarray_i(j, ctrl_pts_opt)[i])
				continue;

			control_pt_index = ctrl_pt_indeces[i];
			get_col_Multiarray_c(j-1, control_points_c)[control_pt_index] += CX_STEP*I;

			update_geo_data_NURBS_parametric_c((const struct const_Multiarray_c*)control_points_c);
			set_up_solver_geometry_c(sim_c);

			dI_dXp_i[col_index] = cimag(functional_c(sim_c))/CX_STEP;

			struct Multiarray_c *R = constructor_residual_c(sensitivity_data, optimization_case);

			for (int k = 0; k < num_res_eqs; k++)
				get_col_Matrix_d(col_index, dR_dXp)[k] = cimag(get_col_Multiarray_c(0, R)[k])/CX_STEP;
			
			destructor_Multiarray_c(R);

			get_col_Multiarray_c(j-1, control_points_c)[control_pt_index] -= CX_STEP*I;

			col_index++;

		}
	}

	//exit(0);

	// Update the geometry one last time with the original control point configurations
	update_geo_data_NURBS_parametric_c((const struct const_Multiarray_c*)control_points_c);
	set_up_solver_geometry_c(sim_c);

	// Free allocated structures
	//destructor_Multiarray_d(R_0);
}
