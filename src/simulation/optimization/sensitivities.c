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

#include "optimization_case.h"
#include "simulation.h"

#include "geometry.h"
#include "geometry_parametric.h"

#include "multiarray.h"
#include "complex_multiarray.h"
#include "multiarray_constructors.h"
#include "complex_multiarray_constructors.h"

#include "solve.h"
#include "solve_implicit.h"

// Static function declarations ************************************************************************************* //

//static struct Multiarray_d* constructor_dR_dXp_finite_diff(struct Optimization_Case *optimization_case);

static void constructor_dI_dXp_and_dR_dXp_finite_diff(struct Optimization_Case *optimization_case);
static struct Multiarray_d* constructor_residual(struct Optimization_Case *optimization_case);

// Interface functions ********************************************************************************************** //

void compute_sensitivities(struct Optimization_Case *optimization_case){

	/*
	Compute the sensitivities of the Residual [dR/dXp] and the objective
	function [dI/dXp] with respect to the design control points Xp. Store the
	results in multiarrays in the optimization_case data structure.

	Arguments:
		optimization_case = Data structure to hold all optimization data

	Return:
		- 
	*/

	constructor_dI_dXp_and_dR_dXp_finite_diff(optimization_case);

}

void destruct_sensitivity_structures(struct Optimization_Case *optimization_case){

	/*
	Free the multiarrays created which hold the sensitivities. This method will
	be called at the end of each design cycle in order to free the dR_dXp and 
	dI_dXp data so as to not have any memory leaks.

	Arguments:
		optimization_case = The optimization_case data structure

	Return:
		-
	*/

	destructor_Multiarray_d(optimization_case->dR_dXp);
	destructor_Multiarray_d(optimization_case->dI_dXp);


}

// Static functions ************************************************************************************************* //

static void constructor_dI_dXp_and_dR_dXp_finite_diff(struct Optimization_Case *optimization_case){

	/*
	Compute the sensitivities of the objective function (I) and residuals (R) with respect to the 
	design variables (Xp). This method computes the components of the dI_dXp structure 
	and dR_dXp structures using finite differences. 

	Ordering of sensitivities
		- If there are n control points on the NURBS patch to use for the optimization, 
			each point may be moved in either the x or y direction (only 2D patches 
			considered for now). Information for the freedom of each point is stored in 
			control_points_optimization in optimization_case->geo_data. 
		- The order is found as follows. For each design point, the partial with respect to the
			x degree of freedom and then the y degree of freedom are listed. This is continued
			for each control point (the ordering of control points is given by the row ordering
			of control_points_optimization)
		- example:
		 	if control_points_optimization = [[1, 0, 1], [2, 1, 1]] then:

			 	The first col of dI_dXp is the partial with respect to Xp1's y degree of freedom, 
			 	the second col is the partial with respect to Xp2's x degree of freedom and 
			 	the third col is the partial with respect to Xp2's y degree of freedom.

	Arguments:
		optimization_case = Data structure to hold all optimization data
	
	*/

	// =================================================
	//                Preprocessing
	// =================================================

	struct Simulation *sim = optimization_case->sim;
	int num_design_pt_dofs = optimization_case->num_design_pts_dofs;
	
	struct Multiarray_i* ctrl_pts_opt = optimization_case->geo_data.control_points_optimization;
	int *ctrl_pt_indeces = get_col_Multiarray_i(0, ctrl_pts_opt);
	int n_pts = (int)ctrl_pts_opt->extents[0];

	// Create the dI_dXp and dR_dXp multiarrays
	struct Multiarray_d* dI_dXp = constructor_empty_Multiarray_d('R',2,(ptrdiff_t[]){1, num_design_pt_dofs});
	double *dI_dXp_i = get_row_Multiarray_d(0, dI_dXp);

	int num_res_eqs = (int)optimization_case->Chi->extents[0];
	struct Multiarray_d* dR_dXp = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_res_eqs, num_design_pt_dofs});
	

	// =================================================
	//     Compute the dI_dXp and dR_dXp components
	// =================================================


	// Compute initial I_0 (objective func) and residual (R_0) values
	double I_0 					= optimization_case->objective_function(sim);
	struct Multiarray_d *R_0 	= constructor_residual(optimization_case);  // free

	struct Multiarray_d* ctrl_pts_and_weights = optimization_case->geo_data.control_points_and_weights;
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
			get_col_Multiarray_d(j-1, ctrl_pts_and_weights)[control_pt_index] += FINITE_DIFF_STEP;


			update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)ctrl_pts_and_weights);
			set_up_solver_geometry(sim);
			dI_dXp_i[col_index] = (optimization_case->objective_function(sim) - I_0)/FINITE_DIFF_STEP;
			
			struct Multiarray_d *R = constructor_residual(optimization_case);  // free
			
			for (int k = 0; k < num_res_eqs; k++)
				get_col_Multiarray_d(col_index, dR_dXp)[k] = 
					(1./FINITE_DIFF_STEP)*(get_col_Multiarray_d(0, R)[k] - get_col_Multiarray_d(0, R_0)[k]);
			
			destructor_Multiarray_d(R);
			

			get_col_Multiarray_d(j-1, ctrl_pts_and_weights)[control_pt_index] -= FINITE_DIFF_STEP;

			col_index++;
			//printf("col_index : %d \n", col_index); fflush(stdout);
		}
	}

	// Update the geometry one last time with the original control point configurations
	update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)ctrl_pts_and_weights);
	set_up_solver_geometry(sim);

	// Free allocated structures
	destructor_Multiarray_d(R_0);


	optimization_case->dI_dXp = dI_dXp;
	optimization_case->dR_dXp = dR_dXp;

}

static struct Multiarray_d* constructor_residual(struct Optimization_Case *optimization_case){

	/*
	Compute the residual vector. Do this by calling compute_rlhs and copying the 
	negative of the RHS vector in the implicit solve.

	Arguments:
		optimization_case = The data structure with all optimization data
	
	Return:
		A multiarray of dimension [num_res_eq x 1] with the value of the resiudal 
		vector.
	*/


	struct Simulation* sim = optimization_case->sim;
	struct Solver_Storage_Implicit* ssi = constructor_Solver_Storage_Implicit(sim); // destructed
	compute_rlhs_optimization(sim, ssi);

	// Set Residual Multiarray
	int num_res_eqs = (int)optimization_case->Chi->extents[0];  // size of adjoint is equivalent to number of residual equations
	struct Multiarray_d* Residual = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_res_eqs, 1});
	double *Residual_i = get_col_Multiarray_d(0, Residual);

	

	// TODO: Find a way to find residual without getting the LHS also (and filling it into
	// PETSc)

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


