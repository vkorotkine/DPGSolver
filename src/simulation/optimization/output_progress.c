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

#include "output_progress.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "macros.h"

#include "simulation.h"
#include "math_functions.h"

#include "matrix.h"
#include "matrix_constructors.h"
#include "matrix_math.h"
#include "multiarray.h"
#include "multiarray_constructors.h"
#include "multiarray_math.h"

#include "objective_functions.h"
#include "optimization_case.h"
#include "gradient.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

FILE* constructor_optimization_progress_file(struct Optimization_Case *optimization_case){

	struct Simulation *sim = optimization_case->sim;

	// Create the optimization convergence output file
	char f_name[4*STRLEN_MAX] = { 0, };
	sprintf(f_name,"%s%c%s%c%s", sim->pde_name,'/',sim->pde_spec,'/',"optimization_convergence.txt");

	char output_name[STRLEN_MAX] = { 0, };
	strcpy(output_name,"../output/");
	strcat(output_name,"paraview/");
	strcat(output_name,f_name);

	FILE* fp;
	if ((fp = fopen(output_name,"w")) == NULL)
		printf("Error: File %s did not open.\n", output_name), exit(1);


	// Print the header for the data
	// - Standard data headers
	fprintf(fp, "Design_Iteration CPU_Time(s) L2_Norm_Gradient_Cost_Function Cost_Function ");

	// - Additional headers (based on the test case)
	if (strstr(optimization_case->optimization_type_spec, "TARGET_CL"))
		fprintf(fp, "CL");

	fprintf(fp, "\n");

	return fp;
}


void destructor_optimization_progress_file(FILE *fp){

	fclose(fp);

}


void progress_file_add_information(FILE *fp, struct Optimization_Case *optimization_case, double L2_grad, 
	double objective_func_value, double time_elapsed, int design_iteration, bool output_to_stdout){


	// \todo
	// In the optimization.data file, have a section where the user can list what data 
	// they want to monitor in the convergence (full names) such as lift coefficient, 
	// pitching moment ... Accordingly, compute those values and output them to the 
	// convergence file.

	// Compute the standard output information and print it

	fprintf(fp, "%d %e %e %e ", design_iteration, time_elapsed, L2_grad, objective_func_value);
	if (output_to_stdout){
		printf("\n************************************************************************************************************\n");
		printf(" Design_Iteration: %d   Time : %f   L2_grad : %e   obj_func : %e   ", 
			design_iteration, time_elapsed, L2_grad, objective_func_value);
	}
	

	// Print additional test specific information
	if (strstr(optimization_case->optimization_type_spec, "TARGET_CL")){

		double CL = compute_Cl(optimization_case->sim);
		
		fprintf(fp, "%e", CL);
			if (output_to_stdout)
				printf("CL: %f ", CL);
	}

	// Print the new line
	fprintf(fp, "\n");
		if (output_to_stdout){
			printf("\n");
			printf("************************************************************************************************************\n\n");
			fflush(stdout);
		}
		
}


void output_NURBS_patch_information(struct Optimization_Case* optimization_case){

	struct Simulation *sim = optimization_case->sim;

	char f_name[4*STRLEN_MAX] = { 0, };
	sprintf(f_name,"%s%c%s%c%s", sim->pde_name,'/',sim->pde_spec,'/',"Optimized_NURBS_Patch.txt");

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

	fprintf(fp, "Control_Point_Data %d \n", (int)optimization_case->geo_data.control_points->extents[0]);
	for (int i = 0; i < (int)optimization_case->geo_data.control_points->extents[0]; i++){
		fprintf(fp, "%.14e %.14e %.14e\n", 
			get_col_Multiarray_d(0, optimization_case->geo_data.control_points)[i],
			get_col_Multiarray_d(1, optimization_case->geo_data.control_points)[i],
			get_col_Multiarray_d(0, optimization_case->geo_data.control_weights)[i]
			);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Control_Point_Connectivity %d %d\n", 
		(int)optimization_case->geo_data.control_pt_wt_connectivity->extents[0],
		(int)optimization_case->geo_data.control_pt_wt_connectivity->extents[1]);
	for (int i = 0; i < (int)optimization_case->geo_data.control_pt_wt_connectivity->extents[0]; i++){
		for (int j = 0; j < (int)optimization_case->geo_data.control_pt_wt_connectivity->extents[1]; j++){
			fprintf(fp, "%d ", get_row_Multiarray_i(i, optimization_case->geo_data.control_pt_wt_connectivity)[j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}


void output_gradient(struct Optimization_Case* optimization_case, 
	double *gradient){

	struct Simulation *sim = optimization_case->sim;
	int n_dof = optimization_case->num_design_pts_dofs; // length of gradient

	char f_name[4*STRLEN_MAX] = { 0, };
	sprintf(f_name,"%s%c%s%c%s", sim->pde_name,'/',sim->pde_spec,'/',"Gradient.txt");

	char output_name[STRLEN_MAX] = { 0, };
	strcpy(output_name,"../output/");
	strcat(output_name,"paraview/");
	strcat(output_name,f_name);

	FILE* fp;
	if ((fp = fopen(output_name,"w")) == NULL)
		printf("Error: File %s did not open.\n", output_name), exit(1);

	for (int i = 0; i < n_dof; i++){
		fprintf(fp, "%f\n", gradient[i]);
	}

	fclose(fp);
}


