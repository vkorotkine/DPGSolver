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

#include "geometry.h"
#include "geometry_parametric.h"


// TODO: Read this from the optimization.data file
#define CONST_GRADIENT_DESCENT_STEP_SIZE 1E-7


// Static function declarations ************************************************************************************* //

static void output_NURBS_patch_information(struct Optimization_Case* optimization_case);

// Interface functions ********************************************************************************************** //


void gradient_descent(struct Optimization_Case *optimization_case, int design_iteration){

	/*
	Traverse in the negative gradient direction to minimize the functional. This function
	will move a small step length in the negative gradient direction. 

	The function will also modify the geometry and call the relevant functions to 
	recompute the metric terms.

	Arguments:
		optimization_case = The optimization case data structure. Holds the information
			about the gradient of the objective.
		design_iteration = The current design iteration.

	Return:
		-
	*/

	UNUSED(design_iteration);

	struct Simulation *sim = optimization_case->sim;

	struct Multiarray_d *ctrl_pts_and_weights = optimization_case->geo_data.control_points_and_weights;
	struct Multiarray_i *ctrl_pts_optimization = optimization_case->geo_data.control_points_optimization;
	struct Multiarray_d *grad_I = optimization_case->grad_I;


	int n_pts_optimization = (int)ctrl_pts_optimization->extents[0];

	int control_pt_index;
	int grad_index = 0;

	for (int i = 0; i < n_pts_optimization; i++){
		// Loop over the design optimization points

		for (int j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point
			// j = 1 (x degree of freeedom) and j = 2 (y degree of freedom)

			if (!get_col_Multiarray_i(j, ctrl_pts_optimization)[i])
				continue;

			// Modify the control point position
			control_pt_index = get_col_Multiarray_i(0, ctrl_pts_optimization)[i];
			get_col_Multiarray_d(j-1, ctrl_pts_and_weights)[control_pt_index] += 
				-1.0 * CONST_GRADIENT_DESCENT_STEP_SIZE * grad_I->data[grad_index++];

		}
	}

	output_NURBS_patch_information(optimization_case);

	// Update the geometry and recompute the metric terms
	update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)ctrl_pts_and_weights);
	set_up_solver_geometry(sim);

	printf("Gradient Descent \n");

}


// Static functions ************************************************************************************************* //


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


