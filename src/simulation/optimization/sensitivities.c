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

#include "macros.h"
#include "definitions_core.h"
#include "definitions_tol.h"

#include "optimization_case.h"
#include "simulation.h"

#include "geometry.h"

#include "multiarray.h"
#include "complex_multiarray.h"
#include "multiarray_constructors.h"
#include "complex_multiarray_constructors.h"

// Static function declarations ************************************************************************************* //

//static struct Multiarray_d* constructor_dR_dXp_finite_diff(struct Optimization_Case *optimization_case);
static struct Multiarray_d* constructor_dI_dXp_finite_diff(struct Optimization_Case *optimization_case);


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

	optimization_case->dI_dXp = constructor_dI_dXp_finite_diff(optimization_case);

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

	UNUSED(optimization_case);
	// TODO: Implement this method

}

// Static functions ************************************************************************************************* //

static struct Multiarray_d* constructor_dI_dXp_finite_diff(struct Optimization_Case *optimization_case){

	/*
	Compute the sensitivities of the objective function (I) with respect to the 
	design variables (Xp). This method computes the components of the dI_dXp structure
	using finite differences.

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

			 	The first term dI_dXp term is the partial with respect to Xp1's y degree of freedom, 
			 	the second term is the partial with respect to Xp2's x degree of freedom and 
			 	the third term is the partial with respect to Xp2's y degree of freedom.

	Arguments:
		optimization_case = Data structure to hold all optimization data
	
	Return:
		Multiarray_d structure with the sensitivities of the objective function
	*/

	// ================================
	//          Preprocessing
	// ================================

	struct Simulation *sim = optimization_case->sim;
	int i,j;

	// - Compute the number of design points dofs (total number of dofs for the
	// 		optimization, which corresponds to the size of the dI_dXp vectors)
	int num_design_pt_dofs = 0;

	struct Multiarray_i* ctrl_pts_opt = optimization_case->geo_data.control_points_optimization;
	int *ctrl_pt_indeces = get_col_Multiarray_i(0, ctrl_pts_opt);
	int n_pts = (int)ctrl_pts_opt->extents[0];
	
	for (i = 0; i < n_pts; i++){
		// Loop over the design points

		for (j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point

			if (get_col_Multiarray_i(j, ctrl_pts_opt)[i])
				num_design_pt_dofs++;
		}
	}

	printf("num_design_pt_dofs : %d \n", num_design_pt_dofs); fflush(stdout);

	struct Multiarray_d* dI_dXp = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_design_pt_dofs,1});
	double *dI_dXp_i = get_col_Multiarray_d(0, dI_dXp);


	// ================================
	//  Compute the dI_dXp components
	// ================================

	// - Perturb each control point and then update the geometry

	int dI_dXp_index = 0;

	double I_0 = optimization_case->objective_function(sim);

	struct Multiarray_d* ctrl_pts_and_weights = optimization_case->geo_data.control_points_and_weights;
	int control_pt_index;

	for (i = 0; i < n_pts; i++){
		// Loop over the design points

		for (j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point
			// j = 1 (x degree of freeedom) and j = 2 (y degree of freedom)

			if (!get_col_Multiarray_i(j, ctrl_pts_opt)[i])
				continue;

			control_pt_index = ctrl_pt_indeces[i];

			get_col_Multiarray_d(j-1, ctrl_pts_and_weights)[control_pt_index] += FINITE_DIFF_STEP;

			// Copy the data into geo_data in geometry_parametric.c
			set_up_solver_geometry(sim);
			dI_dXp_i[dI_dXp_index++] = (optimization_case->objective_function(sim) - I_0)/FINITE_DIFF_STEP;
				
			get_col_Multiarray_d(j-1, ctrl_pts_and_weights)[control_pt_index] -= FINITE_DIFF_STEP;

			printf("dI_dXp_index : %d\n", dI_dXp_index);fflush(stdout);

		}

	}

	return dI_dXp;

}


