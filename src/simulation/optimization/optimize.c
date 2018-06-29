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
#include <time.h>

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
#include "multiarray_constructors.h"
#include "multiarray_math.h"

#include "intrusive.h"
#include "definitions_intrusive.h"

#include "geometry.h"
#include "geometry_parametric.h"
#include "solve_implicit.h"

#include "simulation.h"

#include "volume_solver.h"
#include "face_solver.h"
#include "visualization.h"
#include "definitions_visualization.h"


#include "optimization_case.h"
#include "optimizer_line_search.h"
#include "optimizer_NLPQLP.h"


// TODO: Read from the optimization.data file
#define CONST_L2_GRAD_EXIT 1E-10
#define CONST_OBJECTIVE_FUNC_EXIT 1E-10
#define MAX_NUM_DESIGN_ITERS 250

#define FALSE_FLAG 0

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void optimize(struct Simulation* sim){

	struct Optimization_Case *optimization_case = constructor_Optimization_Case(sim);

	// Optimization routine
	if (strstr(optimization_case->optimizer_spec, "LINE_SEARCH"))
		optimizer_line_search(optimization_case);
	else if(strstr(optimization_case->optimizer_spec, "NLPQLP"))
		optimizer_NLPQLP(optimization_case);
	else
		EXIT_UNSUPPORTED;

	// Clear allocated structures:
	destructor_Optimization_Case(optimization_case);

}




#if FALSE_FLAG

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


// Static functions ************************************************************************************************* //


static void optimize_NLPQLP(struct Optimization_Case *optimization_case){

	/*
	Interface with NLPQLP and use it to perform the shape optimization.

	Arguments:
		optimization_case = The data structure with the optimization data
	*/

	struct Simulation *sim = optimization_case->sim;

	preprocess_optimize_NLPQLP(optimization_case);

	// Optimization routine
	int design_iteration = 0;
	double L2_grad, objective_func_value, CL;

	while (true){

		printf("Design Iteration: %d \n", design_iteration);
	
		// Load the control point locations for the optimizer
		copy_control_pt_data_to_NLPQLP(optimization_case);


		printf("NLPQLP: START, ");

		// Perform the optimization step using NLPQLP
		write_NLPQLP_input_file(optimization_case);  // Write the INPUT.txt file

		// Start a subshell in which the NLPQLP optimizer executable is run
		system("(cd /Users/jm-034232/Documents/McGill/Research/DPGSolver/NLPQLP_Optimizer && ./exec)");

		read_NLPQLP_output_file(optimization_case);  // Process the OUTPUT.txt file

		printf("COMPLETE \n");

		// Process new control point locations from the optimizer
		copy_control_pt_data_from_NLPQLP(optimization_case);  
		


		// Error checking
		if (optimization_case->NLPQLP_data.IFAIL > 0){

			printf("\n ERROR IN OPTIMIZATION : \n");
			printf("IFAIL : %d \n", optimization_case->NLPQLP_data.IFAIL);
			exit(0);
		}


		// Output the optimization progress (for visualization)
		output_NURBS_patch_information(optimization_case);


		// Deform the geometry and solve the flow
		update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)optimization_case->geo_data.control_points_and_weights);
		set_up_solver_geometry(sim);
		solve_implicit(sim);


		// Copy the data into the complex simulation object
		copy_data_r_to_c_sim(optimization_case->sim, optimization_case->sim_c);


		// Break out of the optimization loop (Either the 
		// optimization is succesfully complete or there was an error)
		if (optimization_case->NLPQLP_data.IFAIL == 0)
			break;

		// ============================================
		//    Compute new Objective/Gradient Values
		// ============================================

		// Compute new functional values
		if (optimization_case->NLPQLP_data.IFAIL == -1)
			compute_function_values_NLPQLP(optimization_case);

		// Compute new gradient values
		if (optimization_case->NLPQLP_data.IFAIL == -2){
			compute_gradient_values_NLPQLP(optimization_case);

			// If IFAIL = -2, the design iteration has been completed
			L2_grad = norm_Multiarray_d(optimization_case->NLPQLP_data.dF, "L2");
			objective_func_value = optimization_case->NLPQLP_data.F;
			CL = compute_Cl(sim);

			printf("L2_grad : %e   obj_func : %e   CL : %f \n", L2_grad, objective_func_value, CL);
			
			design_iteration++;

		}
	}

	postprocess_optimize_NLPQLP(optimization_case);

}

static void preprocess_optimize_NLPQLP(struct Optimization_Case *optimization_case){
	
	/*
	Initialize the NLPQLP data and allocate the data structures needed for the NLPQLP
	optimization. Place all data structures in optimization_case.

	Arguments:
		optimization_case = The data structure with the optimization data

	Return:
		-
	*/

	// ========================================
	//       Optimization Parameters
	// ========================================

	// Setup the optimization parameters
	int NP = 1;  // Number of processors

	// Number of design variable dofs
	int N = optimization_case->num_design_pts_dofs;
	int NMAX = N + 1;
	
	// Constrained Optimization Problem with
	// constraint that consecutive design variables must be within delta_y 
	// limit
	
	// NOTE: This only handles the case where all design points move only in the 
	//		y direction (1 dof)

	// Consider an unconstrained optimization problem first

	//int CONST_M = CONST_N-2;
	int M = 0;
	int ME = 0;


	// Independent Optimization Parameters:
	optimization_case->NLPQLP_data.NP 		= NP;
	optimization_case->NLPQLP_data.N 		= N;
	optimization_case->NLPQLP_data.NMAX 	= NMAX;
	optimization_case->NLPQLP_data.M 		= M;
	optimization_case->NLPQLP_data.ME 		= ME;



	// Dependent Optimization Parameters:
	optimization_case->NLPQLP_data.IFAIL 	= 0;  // Initial value is 0
	optimization_case->NLPQLP_data.MODE 	= 0;
	optimization_case->NLPQLP_data.MNN2 	= M + N + N + 2;



	// Remaining Optimization Parameters 
	// Set these up using default values (can be overwritten)
	optimization_case->NLPQLP_data.IOUT 	= 6;
	optimization_case->NLPQLP_data.MAXIT 	= 100;
	optimization_case->NLPQLP_data.MAXFUN 	= 10;
	optimization_case->NLPQLP_data.MAXNM 	= 10;
	optimization_case->NLPQLP_data.LQL 		= 1;
	optimization_case->NLPQLP_data.IPRINT 	= 0;

	optimization_case->NLPQLP_data.ACC 		= 1E-8;
	optimization_case->NLPQLP_data.ACCQP 	= 1E-14;
	optimization_case->NLPQLP_data.STPMIN 	= 1E-10;
	optimization_case->NLPQLP_data.RHO 		= 0.0;



	// Memory Objects (initialize with zeros)
	optimization_case->NLPQLP_data.X 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.X, 0.0);	

	optimization_case->NLPQLP_data.XL 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.XL, 0.0);	

	optimization_case->NLPQLP_data.XU 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.XU, 0.0);	



	optimization_case->NLPQLP_data.F = 0;

	optimization_case->NLPQLP_data.G 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, M});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.G, 0.0);	

	optimization_case->NLPQLP_data.dF 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.dF, 0.0);		

	optimization_case->NLPQLP_data.dG 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){M, NMAX});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.dG, 0.0);	



	optimization_case->NLPQLP_data.C 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){NMAX, NMAX});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.C, 0.0);	

	optimization_case->NLPQLP_data.U 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, optimization_case->NLPQLP_data.MNN2});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.U, 0.0);	

	optimization_case->NLPQLP_data.D 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.D, 0.0);	



	optimization_case->NLPQLP_data.LWA = 3*(optimization_case->NLPQLP_data.N+1)*(optimization_case->NLPQLP_data.N+1)/2 + 
			33*(optimization_case->NLPQLP_data.N+1) + 9*optimization_case->NLPQLP_data.M + 150;
	optimization_case->NLPQLP_data.LWA = optimization_case->NLPQLP_data.LWA*10;	
	
	optimization_case->NLPQLP_data.WA 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, optimization_case->NLPQLP_data.LWA});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.WA, 0.0);	



	optimization_case->NLPQLP_data.LKWA = N + 25;
	optimization_case->NLPQLP_data.LKWA = optimization_case->NLPQLP_data.LKWA*10;	
	
	optimization_case->NLPQLP_data.KWA 	= constructor_empty_Multiarray_i('C',2,(ptrdiff_t[]){1, optimization_case->NLPQLP_data.LKWA});
	set_to_value_Multiarray_i(optimization_case->NLPQLP_data.KWA, 0.0);	



	optimization_case->NLPQLP_data.LACTIV = 2*optimization_case->NLPQLP_data.M + 10;
	optimization_case->NLPQLP_data.ACTIVE   = constructor_empty_Multiarray_i('C',2,(ptrdiff_t[]){1, optimization_case->NLPQLP_data.LACTIV});
	set_to_value_Multiarray_i(optimization_case->NLPQLP_data.ACTIVE, 0.0);	


	// ========================================
	//         Load Optimization Data
	// ========================================

	// Load the design point data
	struct Multiarray_d* ctrl_pts_and_w = optimization_case->geo_data.control_points_and_weights;
	struct Multiarray_i* ctrl_pts_opt = optimization_case->geo_data.control_points_optimization;
	struct Multiarray_d* ctrl_pts_lims = optimization_case->geo_data.control_points_optimization_lims;
	int n_pts = (int)ctrl_pts_opt->extents[0];

	int ctrl_pt_index;
	int vec_index = 0;

	for (int i = 0; i < n_pts; i++){
		// Loop over the design points

		for (int j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point

			if (!get_col_Multiarray_i(j, ctrl_pts_opt)[i])
				continue;

			ctrl_pt_index = get_col_Multiarray_i(0, ctrl_pts_opt)[i];

			optimization_case->NLPQLP_data.X->data[vec_index] = get_col_Multiarray_d(j-1, ctrl_pts_and_w)[ctrl_pt_index];
			optimization_case->NLPQLP_data.XL->data[vec_index] = get_col_Multiarray_d(0, ctrl_pts_lims)[vec_index];
			optimization_case->NLPQLP_data.XU->data[vec_index] = get_col_Multiarray_d(1, ctrl_pts_lims)[vec_index];

			vec_index++;
		}
	}

	// Set the initial Hessian approximation to be a large multiple of the 
	// identity matrix to keep the step length small
	for (int i = 0; i < optimization_case->NLPQLP_data.NMAX; i++){
		get_col_Multiarray_d(i, optimization_case->NLPQLP_data.C)[i] = 500.;
	}
	optimization_case->NLPQLP_data.MODE = 1;


	// Set the initial function values and gradient values
	compute_function_values_NLPQLP(optimization_case);
	compute_gradient_values_NLPQLP(optimization_case);

}


static void postprocess_optimize_NLPQLP(struct Optimization_Case *optimization_case){

	/*
	Postprocess the NLPQLP data and deallocate the data structures created for the NLPQLP
	optimization.

	Arguments:
		optimization_case = The data structure with the optimization data

	Return:
		-
	*/

	destructor_Multiarray_d(optimization_case->NLPQLP_data.X);
	destructor_Multiarray_d(optimization_case->NLPQLP_data.XL);
	destructor_Multiarray_d(optimization_case->NLPQLP_data.XU);
	destructor_Multiarray_d(optimization_case->NLPQLP_data.G);
	destructor_Multiarray_d(optimization_case->NLPQLP_data.U);
	destructor_Multiarray_d(optimization_case->NLPQLP_data.D);
	destructor_Multiarray_d(optimization_case->NLPQLP_data.WA);

	destructor_Multiarray_i(optimization_case->NLPQLP_data.KWA);
	destructor_Multiarray_i(optimization_case->NLPQLP_data.ACTIVE);

}

static void compute_function_values_NLPQLP(struct Optimization_Case *optimization_case){

	/*
	Compute the objective function value and constraint function values. Load the data
	into optimization_case NLPQLP data structure.

	Arguments:
		optimization_case = The data structure with the optimization data

	Return:
		-
	*/

	// Objective Function:
	optimization_case->NLPQLP_data.F = optimization_case->objective_function(optimization_case->sim);

	// Constraint Functions:

	// TODO: Implement constraint functions to be added to the optimization

}

static void compute_gradient_values_NLPQLP(struct Optimization_Case *optimization_case){

	/*
	Compute the gradient of the objective function and constraint functions. Load the data
	into optimization_case NLPQLP data structure.

	Arguments:
		optimization_case = The data structure with the optimization data

	Return:
		-
	*/

	// =================================
	//    Objective Function Gradient
	// =================================

	setup_adjoint(optimization_case);
	solve_adjoint(optimization_case);


	// Compute dI_dXp and dR_dXp
	compute_sensitivities(optimization_case);

	// Use the adjoint and sensitivities to find the gradient of the objective function
	struct Multiarray_d* grad_I = compute_gradient(optimization_case);  // free

	// Load the gradient data into the NLPQLP data structure
	for (int i = 0; i < optimization_case->num_design_pts_dofs; i++)
		optimization_case->NLPQLP_data.dF->data[i] = grad_I->data[i];

	destructor_Multiarray_d(grad_I);
	destruct_sensitivity_structures(optimization_case);
	
	// =================================
	//    Constraint Function Gradient
	// =================================

	// TODO: Implement constraint functions for the optimization

}


static void write_NLPQLP_input_file(struct Optimization_Case *optimization_case){

	/*
	Write the file that will be read by the NLPQLP program.

	NOTE: Uses absolute paths to the input and output files here 
		in the optimization directory
	TODO: Read the path to the NLPQLP directory from an optimization.data file

	Arguments:
		optimization_case = The data structure with the optimization information

	Return:
		-
	*/

	char output_name[STRLEN_MAX] = { 0, };
	strcpy(output_name,"../../NLPQLP_Optimizer/");
	strcat(output_name,"INPUT.txt");

	FILE *fp;

	char lql_string[100];
	int i, j;

	if ((fp = fopen(output_name,"w")) == NULL)
		printf("Error: File %s did not open.\n", output_name), exit(1);

	// Optimizer Properties
	fprintf(fp, "NP \t N \t NMAX \t M \t ME \t IFAIL \t MODE \n");
	fprintf(fp, "%d \t %d \t %d \t %d \t %d \t %d \t %d \n", 
		optimization_case->NLPQLP_data.NP, optimization_case->NLPQLP_data.N, 
		optimization_case->NLPQLP_data.NMAX, optimization_case->NLPQLP_data.M,
		optimization_case->NLPQLP_data.ME, optimization_case->NLPQLP_data.IFAIL, 
		optimization_case->NLPQLP_data.MODE);
	fprintf(fp, "\n");

	fprintf(fp, "LWA \t LKWA \t LACTIV \t IOUT \t ACC \t ACCQP \n");
	fprintf(fp, "%d \t %d \t %d \t %d \t %e \t %e \n", 
		optimization_case->NLPQLP_data.LWA, optimization_case->NLPQLP_data.LKWA, 
		optimization_case->NLPQLP_data.LACTIV, optimization_case->NLPQLP_data.IOUT, 
		optimization_case->NLPQLP_data.ACC, optimization_case->NLPQLP_data.ACCQP);
	fprintf(fp, "\n");

	if(optimization_case->NLPQLP_data.LQL){
		strcpy(lql_string, "T");
	} else{
		strcpy(lql_string, "F");
	}

	fprintf(fp, "STPMIN \t MAXIT \t MAXFUN \t MAXNM \t RHO \t LQL \t IPRINT \n");
	fprintf(fp, "%e \t %d \t %d \t %d \t %e \t %s \t %d \n", 
		optimization_case->NLPQLP_data.STPMIN, optimization_case->NLPQLP_data.MAXIT, 
		optimization_case->NLPQLP_data.MAXFUN, optimization_case->NLPQLP_data.MAXNM,
		optimization_case->NLPQLP_data.RHO, lql_string, optimization_case->NLPQLP_data.IPRINT);
	fprintf(fp, "\n");

	// Design Variable Values
	fprintf(fp, "X \t XL \t XU \n");
	for (i = 0; i < optimization_case->NLPQLP_data.NMAX; i++){
		fprintf(fp, "%.14e \t %.14e \t %.14e \n", 
			optimization_case->NLPQLP_data.X->data[i], 
			optimization_case->NLPQLP_data.XL->data[i], 
			optimization_case->NLPQLP_data.XU->data[i]);
	}
	fprintf(fp, "\n");

	// Objective Function Evaluation
	fprintf(fp, "F \n");
	fprintf(fp, "%.14e \n", optimization_case->NLPQLP_data.F);
	fprintf(fp, "\n");

	// Constraint Function Evaluations
	fprintf(fp, "G \n");
	for (i = 0; i < optimization_case->NLPQLP_data.M; i++)
		fprintf(fp, "%.14e \n", optimization_case->NLPQLP_data.G->data[i]);
	fprintf(fp, "\n");

	// Gradient Objective Function
	fprintf(fp, "dF \n");
	for (i = 0; i < optimization_case->NLPQLP_data.NMAX; i++)
		fprintf(fp, "%.14e ", optimization_case->NLPQLP_data.dF->data[i]);
	fprintf(fp, "\n \n");

	// Gradient Constraint Functions (Column Major Ordering)
	fprintf(fp, "dG \n");
	for (i = 0; i < optimization_case->NLPQLP_data.M; i++){
		for (j = 0; j < optimization_case->NLPQLP_data.NMAX; j++){

			fprintf(fp, "%.14e ", get_col_Multiarray_d(j, optimization_case->NLPQLP_data.dG)[i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	// Hessian Approximation (C)
	fprintf(fp, "C \n");
	for (i = 0; i < optimization_case->NLPQLP_data.NMAX; i++){
		for (j = 0; j < optimization_case->NLPQLP_data.NMAX; j++){
			fprintf(fp, "%.14e ", get_col_Multiarray_d(j, optimization_case->NLPQLP_data.C)[i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	// Multipliers (U)
	fprintf(fp, "U \n");
	for (i = 0; i < optimization_case->NLPQLP_data.MNN2; i++){
		fprintf(fp, "%.14e \n", optimization_case->NLPQLP_data.U->data[i]);
	}
	fprintf(fp, "\n");

	// D Structure
	fprintf(fp, "D \n");
	for (i = 0; i < optimization_case->NLPQLP_data.NMAX; i++){
		fprintf(fp, "%.14e \n", optimization_case->NLPQLP_data.D->data[i]);
	}
	fprintf(fp, "\n");

	// WA Structure
	fprintf(fp, "WA \n");
	for (i = 0; i < optimization_case->NLPQLP_data.LWA; i++){
		fprintf(fp, "%.14e \n", optimization_case->NLPQLP_data.WA->data[i]);
	}
	fprintf(fp, "\n");

	// KWA Structure
	fprintf(fp, "KWA \n");
	for (i = 0; i < optimization_case->NLPQLP_data.LKWA; i++){
		fprintf(fp, "%d \n", optimization_case->NLPQLP_data.KWA->data[i]);
	}
	fprintf(fp, "\n");

	// Active
	fprintf(fp, "ACTIVE \n");
	for (i = 0; i < optimization_case->NLPQLP_data.LACTIV; i++){
		if (optimization_case->NLPQLP_data.ACTIVE->data[i])
			fprintf(fp, "T\n");
		else
			fprintf(fp, "F\n");
	}
	fprintf(fp, "\n");

	fclose(fp);

}


static void read_NLPQLP_output_file(struct Optimization_Case *optimization_case){

	/*
	Read the data from the OUTPUT.txt file (the data that is returned from
	the NLPQLP fortran code). Load the data into the NLPQLP data in the 
	optimization_case data structure.

	Arguments:
		optimization_case = The data structure with the optimization information
	*/

	FILE *fp;

	char line[5000], lql_string[200], *strings, *stringe;
	int i, j;

	double tmpd;
	long tmpl;

	// Open the output file
	char output_name[STRLEN_MAX] = { 0, };
	strcpy(output_name,"../../NLPQLP_Optimizer/");
	strcat(output_name,"OUTPUT.txt");
	if ((fp = fopen(output_name,"r")) == NULL)
		printf("Error: File %s did not open.\n", output_name), exit(1);

	// Properties Block
	fgets(line, sizeof line, fp); // Header
	fgets(line, sizeof line, fp);
	sscanf(line, "%d %d %d %d %d %d %d \n", 
				&optimization_case->NLPQLP_data.NP,
				&optimization_case->NLPQLP_data.N,
				&optimization_case->NLPQLP_data.NMAX,
				&optimization_case->NLPQLP_data.M,
				&optimization_case->NLPQLP_data.ME,
				&optimization_case->NLPQLP_data.IFAIL,
				&optimization_case->NLPQLP_data.MODE);
	fgets(line, sizeof line, fp);  // Space

	fgets(line, sizeof line, fp); // Header
	fgets(line, sizeof line, fp);
	sscanf(line, "%d %d %d %d %lf %lf \n", 
				&optimization_case->NLPQLP_data.LWA,
				&optimization_case->NLPQLP_data.LKWA,
				&optimization_case->NLPQLP_data.LACTIV,
				&optimization_case->NLPQLP_data.IOUT,
				&optimization_case->NLPQLP_data.ACC,
				&optimization_case->NLPQLP_data.ACCQP);
	fgets(line, sizeof line, fp);  // Space

	fgets(line, sizeof line, fp); // Header
	fgets(line, sizeof line, fp);
	sscanf(line, "%lf %d %d %d %lf %s %d \n", 
				&optimization_case->NLPQLP_data.STPMIN,
				&optimization_case->NLPQLP_data.MAXIT,
				&optimization_case->NLPQLP_data.MAXFUN,
				&optimization_case->NLPQLP_data.MAXNM,
				&optimization_case->NLPQLP_data.RHO,
				lql_string,
				&optimization_case->NLPQLP_data.IPRINT);
	fgets(line, sizeof line, fp);  // Space

	if (strstr(lql_string, "T"))
		optimization_case->NLPQLP_data.LQL = 1;
	else
		optimization_case->NLPQLP_data.LQL = 0;

	// Design Points
	fgets(line, sizeof line, fp); // Header

	for (i = 0; i < optimization_case->NLPQLP_data.NMAX; i++){
		fgets(line, sizeof line, fp);

		strings = line;

		// Get the three values from the line
		tmpd = strtod(strings, &stringe);
		strings = stringe;
		optimization_case->NLPQLP_data.X->data[i] = tmpd;

		tmpd = strtod(strings, &stringe);
		strings = stringe;
		optimization_case->NLPQLP_data.XL->data[i] = tmpd;

		tmpd = strtod(strings, &stringe);
		strings = stringe;
		optimization_case->NLPQLP_data.XU->data[i] = tmpd;

	}

	fgets(line, sizeof line, fp);  // Space


	// F (Objective Function Value, outdated so does not need to be saved)
	fgets(line, sizeof line, fp);  // Header
	fgets(line, sizeof line, fp); 
	fgets(line, sizeof line, fp);  // Space


	// G (Constraint Function Values, outdated so does not need to be saved)
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimization_case->NLPQLP_data.M; i++)
		fgets(line, sizeof line, fp); 
	fgets(line, sizeof line, fp);  // Space


	// dF (Gradient of the objective function, outdated so does not need to be saved)
	fgets(line, sizeof line, fp);  // Header
	fgets(line, sizeof line, fp); 
	fgets(line, sizeof line, fp);  // Space


	// dG (Gradient of the constraint functions, outdated so does not need to be saved)
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimization_case->NLPQLP_data.M; i++)
		fgets(line, sizeof line, fp); 
	fgets(line, sizeof line, fp);  // Space


	// C 
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimization_case->NLPQLP_data.NMAX; i++){
		fgets(line, sizeof line, fp);

		strings = line;
		for (j = 0; j < optimization_case->NLPQLP_data.NMAX; j++){
			tmpd = strtod(strings, &stringe);
			strings = stringe;
			get_col_Multiarray_d(j,optimization_case->NLPQLP_data.C)[i] = tmpd;

		}
	}
	fgets(line, sizeof line, fp);  // Space

	// U 
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimization_case->NLPQLP_data.MNN2; i++){
		fgets(line, sizeof line, fp);

		strings = line;
		tmpd = strtod(strings, &stringe);
		strings = stringe;

		optimization_case->NLPQLP_data.U->data[i] = tmpd;

	}
	fgets(line, sizeof line, fp);  // Space


	// D
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimization_case->NLPQLP_data.NMAX; i++){
		fgets(line, sizeof line, fp);

		strings = line;
		tmpd = strtod(strings, &stringe);
		strings = stringe;

		optimization_case->NLPQLP_data.D->data[i] = tmpd;
	}
	fgets(line, sizeof line, fp);  // Space


	// WA
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimization_case->NLPQLP_data.LWA; i++){
		fgets(line, sizeof line, fp);

		strings = line;
		tmpd = strtod(strings, &stringe);
		strings = stringe;

		optimization_case->NLPQLP_data.WA->data[i] = tmpd;
	}
	fgets(line, sizeof line, fp);  // Space


	// KWA
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimization_case->NLPQLP_data.LKWA; i++){
		fgets(line, sizeof line, fp);

		strings = line;
		tmpl = strtol(strings, &stringe, 10);
		strings = stringe;

		optimization_case->NLPQLP_data.KWA->data[i] = (int)tmpl;
	}
	fgets(line, sizeof line, fp);  // Space


	// ACTIVE
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimization_case->NLPQLP_data.LACTIV; i++){
		fgets(line, sizeof line, fp);

		if (strstr(line, "T"))
			optimization_case->NLPQLP_data.ACTIVE->data[i] = 1;
		else
			optimization_case->NLPQLP_data.ACTIVE->data[i] = 0;

	}
	fgets(line, sizeof line, fp);  // Space

	fclose(fp);
}


static void copy_control_pt_data_to_NLPQLP(struct Optimization_Case *optimization_case){

	/*
	This function will be used to transfer the control point data into the NLPQLP
	data structures.

	Arguments:
		optimization_case = The optimization data structure in which data will be 
			modified.
	Return:
		-	

	*/

	struct Multiarray_d* ctrl_pts_and_w = optimization_case->geo_data.control_points_and_weights;
	struct Multiarray_i* ctrl_pts_opt = optimization_case->geo_data.control_points_optimization;
	int n_pts = (int)ctrl_pts_opt->extents[0];

	int ctrl_pt_index;
	int vec_index = 0;

	for (int i = 0; i < n_pts; i++){
		// Loop over the design points

		for (int j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point

			if (!get_col_Multiarray_i(j, ctrl_pts_opt)[i])
				continue;

			// This point can move

			ctrl_pt_index = get_col_Multiarray_i(0, ctrl_pts_opt)[i];

			optimization_case->NLPQLP_data.X->data[vec_index] = get_col_Multiarray_d(j-1, ctrl_pts_and_w)[ctrl_pt_index];
			
			vec_index++;
		}
	}

}


static void copy_control_pt_data_from_NLPQLP(struct Optimization_Case *optimization_case){

	/*
	This function will be used to transfer the control point data from the NLPQLP
	data structures

	TODO: Merge this with the copy to function

	Arguments:
		optimization_case = The optimization data structure in which data will be 
			modified.
	Return:
		-
	*/

	struct Multiarray_d* ctrl_pts_and_w = optimization_case->geo_data.control_points_and_weights;
	struct Multiarray_i* ctrl_pts_opt = optimization_case->geo_data.control_points_optimization;
	int n_pts = (int)ctrl_pts_opt->extents[0];

	int ctrl_pt_index;
	int vec_index = 0;

	for (int i = 0; i < n_pts; i++){
		// Loop over the design control points

		for (int j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point

			if (!get_col_Multiarray_i(j, ctrl_pts_opt)[i])
				continue;

			ctrl_pt_index = get_col_Multiarray_i(0, ctrl_pts_opt)[i];
			
			get_col_Multiarray_d(j-1, ctrl_pts_and_w)[ctrl_pt_index] = optimization_case->NLPQLP_data.X->data[vec_index];

			vec_index++;
		}
	}


}


static void optimize_line_search_method(struct Optimization_Case* optimization_case){

	/*
	Use a line search method (gradient descent or BFGS) to find the optimal shape

	Arguments:
		optimization_case = The data structure holding all optimization data

	Return:
		-
	*/

	struct Simulation *sim = optimization_case->sim;

	preprocessor_minimizer(optimization_case);

	int design_iteration = 0;
	double L2_grad, objective_func_value, CL;

	objective_func_value = optimization_case->objective_function(sim);
	printf(" INITIAL -> obj_func : %e  , Cl: %f \n", objective_func_value, compute_Cl(sim));

	// TODO: Place the file printing in a separate file. Perhaps create an array to 
	// hold up to till max iter values and then just print everything at once

	// Create the optimization convergence output file
	char f_name[4*STRLEN_MAX] = { 0, };
	sprintf(f_name,"%s%c%s%c%s", sim->pde_name,'/',sim->pde_spec,'/',"convergence.txt");

	char output_name[STRLEN_MAX] = { 0, };
	strcpy(output_name,"../output/");
	strcat(output_name,"paraview/");
	strcat(output_name,f_name);

	FILE* fp;
	if ((fp = fopen(output_name,"w")) == NULL)
		printf("Error: File %s did not open.\n", output_name), exit(1);

	// Print the headers for the file
	// Note, we are printing an extra header for CL here. 
	fprintf(fp, "Design_Iteration CPU_Time(s) L2_Norm_Gradient_Cost_Function Cost_Function CL\n");

	// Keep track of the time taken for the optimization
	clock_t start_t, iter_t;
	double t_elapse;
	start_t = clock();

	// ================================
	//      Optimization Routine
	// ================================

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


		// Find the search direction, step length. Then modify the geometry and solve the flow
		//gradient_descent(optimization_case, design_iteration);
		BFGS_minimizer(optimization_case, design_iteration);


		// Copy the new data from the real to the complex structures (for new complex)
		copy_data_r_to_c_sim(optimization_case->sim, optimization_case->sim_c);


		// Compute optimization values to keep track and for the exit condition.
		// Print progress to convergence file
		L2_grad = norm_Multiarray_d(optimization_case->grad_I, "L2");
		objective_func_value = optimization_case->objective_function(sim);
		CL = compute_Cl(sim);  // CL Printing temporary and test specific (TODO fix this)
		iter_t = clock();
		t_elapse = (double)(iter_t - start_t)/(CLOCKS_PER_SEC);

		printf("Time : %f    L2_grad : %e   obj_func : %e   CL : %f \n", t_elapse, L2_grad, objective_func_value, CL);
		fprintf(fp, "%d %e %e %e %e\n", 
			design_iteration+1,
			t_elapse,
			L2_grad,
			objective_func_value,
			CL
			);

		// Destruct allocated data structures:
		destruct_sensitivity_structures(optimization_case);
		destructor_Multiarray_d(optimization_case->grad_I);


		// Exit condition
		if (L2_grad < CONST_L2_GRAD_EXIT || 
			objective_func_value < CONST_OBJECTIVE_FUNC_EXIT ||
			design_iteration >= MAX_NUM_DESIGN_ITERS)
			break;

		design_iteration++;

	}

	// ================================
	//         Postprocessing
	// ================================

	postprocessor_minimizer(optimization_case);
	fclose(fp);

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

	int num_design_dofs = optimization_case->num_design_pts_dofs;

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

#endif



