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

#include "optimizer_NLPQLP.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#include "macros.h"
#include "simulation.h"
#include "file_processing.h"
#include "definitions_alloc.h"

#include "matrix.h"
#include "matrix_constructors.h"
#include "matrix_math.h"
#include "multiarray.h"
#include "multiarray_constructors.h"
#include "multiarray_math.h"

#include "geometry.h"
#include "geometry_parametric.h"
#include "solve_implicit.h"

#include "math_functions.h"

#include "optimization_case.h"
#include "adjoint.h"
#include "sensitivities.h"
#include "gradient.h"
#include "output_progress.h"

// TESTING
#include "functionals.h"

// Static function declarations ************************************************************************************* //


/** \brief Construct the data structure to hold the NLPQLP information
 */
static struct Optimizer_NLPQLP_Data* constructor_Optimizer_NLPQLP_Data(
	struct Optimization_Case* optimization_case ///< Consult optimization_case.h
	);


/** \brief Destruct the data structure that holds the optimization NLPQLP information along
 * with the allocated memory that it holds.
 */
static void destructor_Optimizer_NLPQLP_Data(
	struct Optimizer_NLPQLP_Data* optimizer_nlpqlp_data ///< Consult optimizer_NLPQLP.h
	);


/** \brief 	Read the required data from the optimization.geo file and load it into the 
 * 	Optimizer_NLPQLP_Data data structure.
 */
static void read_optimization_data(
	struct Optimizer_NLPQLP_Data* optimizer_nlpqlp_data ///< Consult optimizer_NLPQLP.h
	);


/** \brief Compute the objective function value and constraint function values (if applicable). 
 * 	Load the data into the Optimizer_NLPQLP_Data data structure. 
 */
static void compute_function_values_NLPQLP(
	struct Optimization_Case* optimization_case, ///< Consult optimization_case.h
	struct Optimizer_NLPQLP_Data *optimizer_nlpqlp_data ///< Consult optimizer_NLPQLP.h
	);


/** \brief 	Compute the gradient of the objective function and constraint functions. Load the data
 *	into Optimizer_NLPQLP_Data data structure.
 */
static void compute_gradient_values_NLPQLP(
	struct Optimization_Case *optimization_case, ///< Consult optimization_case.h
	struct Optimizer_NLPQLP_Data* optimizer_nlpqlp_data ///< Consult optimizer_NLPQLP.h
	);


/** \brief Write the file that will be read by the NLPQLP program.
 *
 *	NOTE: Uses absolute paths to the input and output files here 
 *		in the optimization directory
 */
static void write_NLPQLP_input_file(
	struct Optimizer_NLPQLP_Data* optimizer_nlpqlp_data ///< Consult optimizer_NLPQLP.h
	);


/** \brief Read the data from the OUTPUT.txt file (the data that is returned from
 *	the NLPQLP fortran executable). Load the data into the NLPQLP data structure
 */
static void read_NLPQLP_output_file(
	struct Optimizer_NLPQLP_Data* optimizer_nlpqlp_data ///< Consult optimizer_NLPQLP.h
	);


/** \brief Function used to transfer the control point data to and from the NLPQLP
 *	data structures.
 */
static void transfer_control_point_data_NLPQLP(
	struct Optimization_Case *optimization_case, ///< Consult optimization_case.h
	struct Optimizer_NLPQLP_Data* optimizer_nlpqlp_data, ///< Consult optimizer_NLPQLP.h
	char transfer_direction ///< 'o' = transfer to 'o'ptimiziation_case data structure, 'n' = transfer to 'n'lpqlp data structure
	);

// Interface functions ********************************************************************************************** //

void optimizer_NLPQLP(struct Optimization_Case* optimization_case){

	// ================================
	//         Preprocessing
	// ================================

	struct Simulation *sim = optimization_case->sim;

	FILE *fp = constructor_optimization_progress_file(optimization_case);
	struct Optimizer_NLPQLP_Data *optimizer_nlpqlp_data = constructor_Optimizer_NLPQLP_Data(optimization_case);

	// Keep track of the time taken for the optimization and number of design iterations
	clock_t start_t, iter_t;
	double t_elapse;
	start_t = clock();

	int design_iteration = 0;
	bool print_initial_values = true;

	// ================================
	//      Optimization Routine
	// ================================

	while (true){
	
		// Print the progress for the very first design iteration (the initial values before
		// optimization has begun).
		if (print_initial_values){
			print_initial_values = false;

			double L2_grad = norm_Multiarray_d(optimizer_nlpqlp_data->dF, "L2");
			double objective_func_value = optimization_case->objective_function(sim);
			iter_t = clock();
			t_elapse = (double)(iter_t - start_t)/(CLOCKS_PER_SEC);

			progress_file_add_information(fp, optimization_case, L2_grad, objective_func_value, 
				t_elapse, design_iteration, true);

			// Output the initial gradient
			output_gradient(optimization_case, optimizer_nlpqlp_data->dF->data);
		}

		printf("\nStart NLPQLP  ");

		char transfer_key;

		// Load the control point locations for the optimizer
		transfer_key = 'n';
		transfer_control_point_data_NLPQLP(optimization_case, optimizer_nlpqlp_data, transfer_key);

		// Perform the optimization step using NLPQLP
		write_NLPQLP_input_file(optimizer_nlpqlp_data);  // Write the INPUT.txt file

		// Start a subshell in which the NLPQLP optimizer executable will be run.		
		// Set the executaion command first and then start the subshell
		char nlpqlp_execution_command[STRLEN_MAX*4];
		int index = 0;
		index += sprintf(nlpqlp_execution_command, "(cd ");
		index += sprintf(nlpqlp_execution_command + index, "%s", optimizer_nlpqlp_data->nlpqlp_directory_abs_path);
		index += sprintf(nlpqlp_execution_command + index, " && ./exec)");
		system(nlpqlp_execution_command);

		read_NLPQLP_output_file(optimizer_nlpqlp_data);  // Process the OUTPUT.txt file

		// Process new control point locations from the optimizer
		transfer_key = 'o';
		transfer_control_point_data_NLPQLP(optimization_case, optimizer_nlpqlp_data, transfer_key);
		
		printf("Completed NLPQLP : IFAIL = %d \n\n", optimizer_nlpqlp_data->IFAIL);

		// Error checking
		if (optimizer_nlpqlp_data->IFAIL > 0){

			printf("\n ERROR IN OPTIMIZATION : \n");
			printf("IFAIL : %d \n", optimizer_nlpqlp_data->IFAIL);
			break; // break because still want to print data to file
		}


		// Output the optimization progress (for visualization)
		output_NURBS_patch_information(optimization_case);


		// Deform the geometry and solve the flow
		update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)optimization_case->geo_data.control_points, (const struct Simulation *) sim);
		set_up_solver_geometry(sim);
		solve_implicit(sim);


		// With the updated geometry and flow, copy the data into the complex simulation object
		copy_data_r_to_c_sim(optimization_case->sim, optimization_case->sim_c);
		destructor_Multiarray_c(optimization_case->geo_data.control_points_c);
		optimization_case->geo_data.control_points_c = 
			constructor_copy_Multiarray_c_Multiarray_d(optimization_case->geo_data.control_points);


		// Optimization successful so break out of the loop
		if (optimizer_nlpqlp_data->IFAIL == 0)
			break;


		// ============================================
		//    Compute new Functional/Gradient Values
		// ============================================


		// Compute new functional values
		if (optimizer_nlpqlp_data->IFAIL == -1)
			compute_function_values_NLPQLP(optimization_case, optimizer_nlpqlp_data);


		// Compute new gradient values
		if (optimizer_nlpqlp_data->IFAIL == -2){
			
			// If IFAIL = -2, the design iteration has been completed
			design_iteration++;

			compute_gradient_values_NLPQLP(optimization_case, optimizer_nlpqlp_data);

			// Output the progress for the recently completed design iteration
			double L2_grad = norm_Multiarray_d(optimizer_nlpqlp_data->dF, "L2");
			double objective_func_value = optimization_case->objective_function(sim);
			iter_t = clock();
			t_elapse = (double)(iter_t - start_t)/(CLOCKS_PER_SEC);
			
			progress_file_add_information(fp, optimization_case, L2_grad, objective_func_value, 
				t_elapse, design_iteration, true);

		}
	}

	// ================================
	//         Postprocessing
	// ================================

	destructor_optimization_progress_file(fp);
	destructor_Optimizer_NLPQLP_Data(optimizer_nlpqlp_data);
}


// Static functions ************************************************************************************************* //


static struct Optimizer_NLPQLP_Data* constructor_Optimizer_NLPQLP_Data(struct Optimization_Case* optimization_case){

	struct Optimizer_NLPQLP_Data* optimizer_nlpqlp_data = calloc(1,sizeof *optimizer_nlpqlp_data); // returned

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

	int M = optimization_case->num_total_constraints;
	int ME = optimization_case->num_equality_constraints;


	// Independent Optimization Parameters:
	optimizer_nlpqlp_data->NP 		= NP;
	optimizer_nlpqlp_data->N 		= N;
	optimizer_nlpqlp_data->NMAX 	= NMAX;
	optimizer_nlpqlp_data->M 		= M;
	optimizer_nlpqlp_data->ME 		= ME;



	// Dependent Optimization Parameters:
	optimizer_nlpqlp_data->IFAIL 	= 0;  // Initial value is 0
	optimizer_nlpqlp_data->MODE 	= 0;
	optimizer_nlpqlp_data->MNN2 	= M + N + N + 2;


	// Settable Optimization Parameters 
	read_optimization_data(optimizer_nlpqlp_data);


	// Memory Objects (initialize with zeros)
	optimizer_nlpqlp_data->X 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimizer_nlpqlp_data->X, 0.0);	

	optimizer_nlpqlp_data->XL 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimizer_nlpqlp_data->XL, 0.0);	

	optimizer_nlpqlp_data->XU 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimizer_nlpqlp_data->XU, 0.0);	



	optimizer_nlpqlp_data->F = 0;

	optimizer_nlpqlp_data->G 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, M});
	set_to_value_Multiarray_d(optimizer_nlpqlp_data->G, 0.0);	

	optimizer_nlpqlp_data->dF 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimizer_nlpqlp_data->dF, 0.0);		

	optimizer_nlpqlp_data->dG 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){M, NMAX});
	set_to_value_Multiarray_d(optimizer_nlpqlp_data->dG, 0.0);	



	optimizer_nlpqlp_data->C 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){NMAX, NMAX});
	set_to_value_Multiarray_d(optimizer_nlpqlp_data->C, 0.0);	

	optimizer_nlpqlp_data->U 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, optimizer_nlpqlp_data->MNN2});
	set_to_value_Multiarray_d(optimizer_nlpqlp_data->U, 0.0);	

	optimizer_nlpqlp_data->D 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimizer_nlpqlp_data->D, 0.0);	



	optimizer_nlpqlp_data->LWA = 3*(optimizer_nlpqlp_data->N+1)*(optimizer_nlpqlp_data->N+1)/2 + 
			33*(optimizer_nlpqlp_data->N+1) + 9*optimizer_nlpqlp_data->M + 150;
	optimizer_nlpqlp_data->LWA = optimizer_nlpqlp_data->LWA*10;	
	
	optimizer_nlpqlp_data->WA 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, optimizer_nlpqlp_data->LWA});
	set_to_value_Multiarray_d(optimizer_nlpqlp_data->WA, 0.0);	



	optimizer_nlpqlp_data->LKWA = N + 25;
	optimizer_nlpqlp_data->LKWA = optimizer_nlpqlp_data->LKWA*10;	
	
	optimizer_nlpqlp_data->KWA 	= constructor_empty_Multiarray_i('C',2,(ptrdiff_t[]){1, optimizer_nlpqlp_data->LKWA});
	set_to_value_Multiarray_i(optimizer_nlpqlp_data->KWA, 0.0);	



	optimizer_nlpqlp_data->LACTIV = 2*optimizer_nlpqlp_data->M + 10;
	optimizer_nlpqlp_data->ACTIVE   = constructor_empty_Multiarray_i('C',2,(ptrdiff_t[]){1, optimizer_nlpqlp_data->LACTIV});
	set_to_value_Multiarray_i(optimizer_nlpqlp_data->ACTIVE, 0.0);	


	// ========================================
	//         Load Optimization Data
	// ========================================

	// Load the design point data
	// From here thru end of for loop def needs to change to accomodate multipatch
	if(optimization_case->sim->nurbs_multipatch){
		for(int patch_index=0;patch_index<optimization_case->num_patches;patch_index++){
			struct Multiarray_d* ctrl_pts = optimization_case->NURBS_Patch_Data[patch_index].control_points;
			struct Multiarray_i* ctrl_pts_opt = optimization_case->NURBS_Patch_Data[patch_index].control_points_optimization;
			struct Multiarray_d* ctrl_pts_lims = optimization_case->NURBS_Patch_Data[patch_index].control_points_optimization_lims;
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

					optimizer_nlpqlp_data->X->data[vec_index]  = get_col_Multiarray_d(j-1, ctrl_pts)[ctrl_pt_index];
					optimizer_nlpqlp_data->XL->data[vec_index] = get_col_Multiarray_d(0, ctrl_pts_lims)[vec_index];
					optimizer_nlpqlp_data->XU->data[vec_index] = get_col_Multiarray_d(1, ctrl_pts_lims)[vec_index];

					vec_index++;
				}
			}
		}
	}
	else{
		struct Multiarray_d* ctrl_pts = optimization_case->geo_data.control_points;
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

				optimizer_nlpqlp_data->X->data[vec_index]  = get_col_Multiarray_d(j-1, ctrl_pts)[ctrl_pt_index];
				optimizer_nlpqlp_data->XL->data[vec_index] = get_col_Multiarray_d(0, ctrl_pts_lims)[vec_index];
				optimizer_nlpqlp_data->XU->data[vec_index] = get_col_Multiarray_d(1, ctrl_pts_lims)[vec_index];

				vec_index++;
			}
		}
	}

	// Set the initial Hessian approximation to be a large multiple of the 
	// identity matrix to keep the step length small
	for (int i = 0; i < optimizer_nlpqlp_data->NMAX; i++){
		get_col_Multiarray_d(i, optimizer_nlpqlp_data->C)[i] = 500.;
	}
	optimizer_nlpqlp_data->MODE = 1;

	// Set the initial function values and gradient values
	compute_function_values_NLPQLP(optimization_case, optimizer_nlpqlp_data);
	compute_gradient_values_NLPQLP(optimization_case, optimizer_nlpqlp_data);

	return optimizer_nlpqlp_data;
}


static void destructor_Optimizer_NLPQLP_Data(struct Optimizer_NLPQLP_Data* optimizer_nlpqlp_data){

	destructor_Multiarray_d(optimizer_nlpqlp_data->X);
	destructor_Multiarray_d(optimizer_nlpqlp_data->XL);
	destructor_Multiarray_d(optimizer_nlpqlp_data->XU);

	destructor_Multiarray_d(optimizer_nlpqlp_data->G);
	destructor_Multiarray_d(optimizer_nlpqlp_data->dF);
	destructor_Multiarray_d(optimizer_nlpqlp_data->dG);
	
	destructor_Multiarray_d(optimizer_nlpqlp_data->C);
	destructor_Multiarray_d(optimizer_nlpqlp_data->U);
	destructor_Multiarray_d(optimizer_nlpqlp_data->D);
	destructor_Multiarray_d(optimizer_nlpqlp_data->WA);

	destructor_Multiarray_i(optimizer_nlpqlp_data->KWA);
	destructor_Multiarray_i(optimizer_nlpqlp_data->ACTIVE);

	free((void*)optimizer_nlpqlp_data);
}


static void read_optimization_data(struct Optimizer_NLPQLP_Data* optimizer_nlpqlp_data){

	// Get the file pointer to the optimization file
	FILE* input_file = fopen_input('o',NULL,NULL); // closed
	char line[STRLEN_MAX];

	// Read in the information from the file
	while (fgets(line,sizeof(line),input_file)) {

		// Any line with a comment flag should be skipped
		if (strstr(line, "//"))
			continue;

		if (strstr(line, "IOUT")) 	read_skip_i_1(line, 1, &optimizer_nlpqlp_data->IOUT, 1);
		if (strstr(line, "MAXIT")) 	read_skip_i_1(line, 1, &optimizer_nlpqlp_data->MAXIT, 1);
		if (strstr(line, "MAXFUN")) read_skip_i_1(line, 1, &optimizer_nlpqlp_data->MAXFUN, 1);
		if (strstr(line, "MAXNM")) 	read_skip_i_1(line, 1, &optimizer_nlpqlp_data->MAXNM, 1);
		if (strstr(line, "LQL")) 	read_skip_i_1(line, 1, &optimizer_nlpqlp_data->LQL, 1);
		if (strstr(line, "IPRINT")) read_skip_i_1(line, 1, &optimizer_nlpqlp_data->IPRINT, 1);

		if (strstr(line, "ACC")) 	read_skip_d_1(line, 1, &optimizer_nlpqlp_data->ACC, 1);
		if (strstr(line, "ACCQP")) 	read_skip_d_1(line, 1, &optimizer_nlpqlp_data->ACCQP, 1);
		if (strstr(line, "STPMIN")) read_skip_d_1(line, 1, &optimizer_nlpqlp_data->STPMIN, 1);
		if (strstr(line, "RHO")) 	read_skip_d_1(line, 1, &optimizer_nlpqlp_data->RHO, 1);

		if (strstr(line, "optimizer_nlpqlp_abs_path")) read_skip_c(line, optimizer_nlpqlp_data->nlpqlp_directory_abs_path);
	}

	fclose(input_file);
}


static void compute_function_values_NLPQLP(struct Optimization_Case *optimization_case, 
	struct Optimizer_NLPQLP_Data *optimizer_nlpqlp_data){

	struct Simulation *sim = optimization_case->sim;

	// ==================================
	//         Objective Function
	// ==================================

	optimizer_nlpqlp_data->F = optimization_case->objective_function(sim);


	// ==================================
	//        Constraint Functions
	// ==================================

	struct Constraint_Function_Data* constraint_function_data = optimization_case->constraint_function_data;
	int constraint_function_index = 0;
	double constraint_function_value;

	while(true){

		if(constraint_function_data == NULL)
			break;

		// Compute the constraint function value and store it in the G structure
		constraint_function_value =  constraint_function_data->functional_f(sim);
		constraint_function_value += constraint_function_data->k;
		constraint_function_value *= constraint_function_data->a;

		get_col_Multiarray_d(constraint_function_index, optimizer_nlpqlp_data->G)[0] = constraint_function_value;

		constraint_function_data = constraint_function_data->next;		

		constraint_function_index++;
	}
}


static void compute_gradient_values_NLPQLP(struct Optimization_Case *optimization_case, 
	struct Optimizer_NLPQLP_Data* optimizer_nlpqlp_data){

	struct Simulation *sim   = optimization_case->sim;
	struct Simulation *sim_c = optimization_case->sim_c;

	// ===================================
	//     Objective Function Gradient
	// ===================================

	struct Adjoint_Data *adjoint_data 		= constructor_Adjoint_Data(optimization_case);
	struct Sensitivity_Data *sensivity_data = constructor_Sensitivity_Data(optimization_case);
	struct Gradient_Data *gradient_data 	= constructor_Gradient_Data(adjoint_data, sensivity_data);


	// Solve the Adjoint equation
	setup_adjoint(adjoint_data, sim, sim_c, optimization_case->objective_function, optimization_case->objective_function_c);
	solve_adjoint(adjoint_data, sim);


	// Compute the Sensitivities
	compute_sensitivities(sensivity_data, optimization_case, optimization_case->objective_function, optimization_case->objective_function_c);


	// Compute the gradient using the sensitivities and the adjoint
	compute_gradient(gradient_data, adjoint_data, sensivity_data);

	// Store gradient data in NLPQLP data structure
	for (int i = 0; i < optimization_case->num_design_pts_dofs; i++)
		optimizer_nlpqlp_data->dF->data[i] = gradient_data->Gradient->data[i];


	// Destruct allocated data structures:
	destructor_Adjoint_Data(adjoint_data);
	destructor_Sensitivity_Data(sensivity_data);
	destructor_Gradient_Data(gradient_data);


	// ===================================
	//    Constraint Function Gradients
	// ===================================	

	struct Constraint_Function_Data* constraint_function_data = optimization_case->constraint_function_data;
	int constraint_function_index = 0;

	while(true){

		if(constraint_function_data == NULL)
			break;


		// Compute the gradient of the given constraint function by using the adjoint
		// approach.

		struct Adjoint_Data *adjoint_data 		= constructor_Adjoint_Data(optimization_case);
		struct Sensitivity_Data *sensivity_data = constructor_Sensitivity_Data(optimization_case);
		struct Gradient_Data *gradient_data 	= constructor_Gradient_Data(adjoint_data, sensivity_data);


		// Solve the Adjoint equation
		setup_adjoint(adjoint_data, sim, sim_c,	constraint_function_data->functional_f, constraint_function_data->functional_f_c);
		solve_adjoint(adjoint_data, sim);


		// Compute the Sensitivities
		compute_sensitivities(sensivity_data, optimization_case, constraint_function_data->functional_f, constraint_function_data->functional_f_c);


		// Compute the gradient using the sensitivities and the adjoint. Scale the gradients
		// by the specified multiplier (due to the chain rule)
		compute_gradient(gradient_data, adjoint_data, sensivity_data);
		for (int i = 0; i < optimization_case->num_design_pts_dofs; i++)
			gradient_data->Gradient->data[i] *= constraint_function_data->a;


		// Store gradient data in NLPQLP data structure (dG)
		for (int j = 0; j < optimization_case->num_design_pts_dofs; j++){
			get_col_Multiarray_d(j, optimizer_nlpqlp_data->dG)[constraint_function_index] = gradient_data->Gradient->data[j];
		}


		// Destruct allocated data structures:
		destructor_Adjoint_Data(adjoint_data);
		destructor_Sensitivity_Data(sensivity_data);
		destructor_Gradient_Data(gradient_data);


		constraint_function_data = constraint_function_data->next;		
		constraint_function_index++;
	}
}


static void write_NLPQLP_input_file(struct Optimizer_NLPQLP_Data* optimizer_nlpqlp_data){

	char output_name[STRLEN_MAX] = { 0, };
	strcpy(output_name, optimizer_nlpqlp_data->nlpqlp_directory_abs_path);
	strcat(output_name,"/INPUT.txt");

	FILE *fp;

	char lql_string[100];
	int i, j;

	if ((fp = fopen(output_name,"w")) == NULL)
		printf("Error: File %s did not open.\n", output_name), exit(1);

	// Optimizer Properties
	fprintf(fp, "NP \t N \t NMAX \t M \t ME \t IFAIL \t MODE \n");
	fprintf(fp, "%d \t %d \t %d \t %d \t %d \t %d \t %d \n", 
		optimizer_nlpqlp_data->NP, optimizer_nlpqlp_data->N, 
		optimizer_nlpqlp_data->NMAX, optimizer_nlpqlp_data->M,
		optimizer_nlpqlp_data->ME, optimizer_nlpqlp_data->IFAIL, 
		optimizer_nlpqlp_data->MODE);
	fprintf(fp, "\n");

	fprintf(fp, "LWA \t LKWA \t LACTIV \t IOUT \t ACC \t ACCQP \n");
	fprintf(fp, "%d \t %d \t %d \t %d \t %e \t %e \n", 
		optimizer_nlpqlp_data->LWA, optimizer_nlpqlp_data->LKWA, 
		optimizer_nlpqlp_data->LACTIV, optimizer_nlpqlp_data->IOUT, 
		optimizer_nlpqlp_data->ACC, optimizer_nlpqlp_data->ACCQP);
	fprintf(fp, "\n");

	if(optimizer_nlpqlp_data->LQL){
		strcpy(lql_string, "T");
	} else{
		strcpy(lql_string, "F");
	}

	fprintf(fp, "STPMIN \t MAXIT \t MAXFUN \t MAXNM \t RHO \t LQL \t IPRINT \n");
	fprintf(fp, "%e \t %d \t %d \t %d \t %e \t %s \t %d \n", 
		optimizer_nlpqlp_data->STPMIN, optimizer_nlpqlp_data->MAXIT, 
		optimizer_nlpqlp_data->MAXFUN, optimizer_nlpqlp_data->MAXNM,
		optimizer_nlpqlp_data->RHO, lql_string, optimizer_nlpqlp_data->IPRINT);
	fprintf(fp, "\n");

	// Design Variable Values
	fprintf(fp, "X \t XL \t XU \n");
	for (i = 0; i < optimizer_nlpqlp_data->NMAX; i++){
		fprintf(fp, "%.14e \t %.14e \t %.14e \n", 
			optimizer_nlpqlp_data->X->data[i], 
			optimizer_nlpqlp_data->XL->data[i], 
			optimizer_nlpqlp_data->XU->data[i]);
	}
	fprintf(fp, "\n");

	// Objective Function Evaluation
	fprintf(fp, "F \n");
	fprintf(fp, "%.14e \n", optimizer_nlpqlp_data->F);
	fprintf(fp, "\n");

	// Constraint Function Evaluations
	fprintf(fp, "G \n");
	for (i = 0; i < optimizer_nlpqlp_data->M; i++)
		fprintf(fp, "%.14e \n", optimizer_nlpqlp_data->G->data[i]);
	fprintf(fp, "\n");

	// Gradient Objective Function
	fprintf(fp, "dF \n");
	for (i = 0; i < optimizer_nlpqlp_data->NMAX; i++)
		fprintf(fp, "%.14e ", optimizer_nlpqlp_data->dF->data[i]);
	fprintf(fp, "\n \n");

	// Gradient Constraint Functions (Column Major Ordering)
	fprintf(fp, "dG \n");
	for (i = 0; i < optimizer_nlpqlp_data->M; i++){
		for (j = 0; j < optimizer_nlpqlp_data->NMAX; j++){

			fprintf(fp, "%.14e ", get_col_Multiarray_d(j, optimizer_nlpqlp_data->dG)[i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	// Hessian Approximation (C)
	fprintf(fp, "C \n");
	for (i = 0; i < optimizer_nlpqlp_data->NMAX; i++){
		for (j = 0; j < optimizer_nlpqlp_data->NMAX; j++){
			fprintf(fp, "%.14e ", get_col_Multiarray_d(j, optimizer_nlpqlp_data->C)[i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	// Multipliers (U)
	fprintf(fp, "U \n");
	for (i = 0; i < optimizer_nlpqlp_data->MNN2; i++){
		fprintf(fp, "%.14e \n", optimizer_nlpqlp_data->U->data[i]);
	}
	fprintf(fp, "\n");

	// D Structure
	fprintf(fp, "D \n");
	for (i = 0; i < optimizer_nlpqlp_data->NMAX; i++){
		fprintf(fp, "%.14e \n", optimizer_nlpqlp_data->D->data[i]);
	}
	fprintf(fp, "\n");

	// WA Structure
	fprintf(fp, "WA \n");
	for (i = 0; i < optimizer_nlpqlp_data->LWA; i++){
		fprintf(fp, "%.14e \n", optimizer_nlpqlp_data->WA->data[i]);
	}
	fprintf(fp, "\n");

	// KWA Structure
	fprintf(fp, "KWA \n");
	for (i = 0; i < optimizer_nlpqlp_data->LKWA; i++){
		fprintf(fp, "%d \n", optimizer_nlpqlp_data->KWA->data[i]);
	}
	fprintf(fp, "\n");

	// Active
	fprintf(fp, "ACTIVE \n");
	for (i = 0; i < optimizer_nlpqlp_data->LACTIV; i++){
		if (optimizer_nlpqlp_data->ACTIVE->data[i])
			fprintf(fp, "T\n");
		else
			fprintf(fp, "F\n");
	}
	fprintf(fp, "\n");

	fclose(fp);
}


static void read_NLPQLP_output_file(struct Optimizer_NLPQLP_Data* optimizer_nlpqlp_data){

	FILE *fp;

	char line[5000], lql_string[200], *strings, *stringe;
	int i, j;

	double tmpd;
	long tmpl;

	// Open the output file
	char output_name[STRLEN_MAX] = { 0, };
	strcpy(output_name, optimizer_nlpqlp_data->nlpqlp_directory_abs_path);
	strcat(output_name,"/OUTPUT.txt");
	if ((fp = fopen(output_name,"r")) == NULL)
		printf("Error: File %s did not open.\n", output_name), exit(1);

	// Properties Block
	fgets(line, sizeof line, fp); // Header
	fgets(line, sizeof line, fp);
	sscanf(line, "%d %d %d %d %d %d %d \n", 
				&optimizer_nlpqlp_data->NP,
				&optimizer_nlpqlp_data->N,
				&optimizer_nlpqlp_data->NMAX,
				&optimizer_nlpqlp_data->M,
				&optimizer_nlpqlp_data->ME,
				&optimizer_nlpqlp_data->IFAIL,
				&optimizer_nlpqlp_data->MODE);
	fgets(line, sizeof line, fp);  // Space

	fgets(line, sizeof line, fp); // Header
	fgets(line, sizeof line, fp);
	sscanf(line, "%d %d %d %d %lf %lf \n", 
				&optimizer_nlpqlp_data->LWA,
				&optimizer_nlpqlp_data->LKWA,
				&optimizer_nlpqlp_data->LACTIV,
				&optimizer_nlpqlp_data->IOUT,
				&optimizer_nlpqlp_data->ACC,
				&optimizer_nlpqlp_data->ACCQP);
	fgets(line, sizeof line, fp);  // Space

	fgets(line, sizeof line, fp); // Header
	fgets(line, sizeof line, fp);
	sscanf(line, "%lf %d %d %d %lf %s %d \n", 
				&optimizer_nlpqlp_data->STPMIN,
				&optimizer_nlpqlp_data->MAXIT,
				&optimizer_nlpqlp_data->MAXFUN,
				&optimizer_nlpqlp_data->MAXNM,
				&optimizer_nlpqlp_data->RHO,
				lql_string,
				&optimizer_nlpqlp_data->IPRINT);
	fgets(line, sizeof line, fp);  // Space

	if (strstr(lql_string, "T"))
		optimizer_nlpqlp_data->LQL = 1;
	else
		optimizer_nlpqlp_data->LQL = 0;

	// Design Points
	fgets(line, sizeof line, fp); // Header

	for (i = 0; i < optimizer_nlpqlp_data->NMAX; i++){
		fgets(line, sizeof line, fp);

		strings = line;

		// Get the three values from the line
		tmpd = strtod(strings, &stringe);
		strings = stringe;
		optimizer_nlpqlp_data->X->data[i] = tmpd;

		tmpd = strtod(strings, &stringe);
		strings = stringe;
		optimizer_nlpqlp_data->XL->data[i] = tmpd;

		tmpd = strtod(strings, &stringe);
		strings = stringe;
		optimizer_nlpqlp_data->XU->data[i] = tmpd;

	}

	fgets(line, sizeof line, fp);  // Space


	// F (Objective Function Value, outdated so does not need to be saved)
	fgets(line, sizeof line, fp);  // Header
	fgets(line, sizeof line, fp); 
	fgets(line, sizeof line, fp);  // Space


	// G (Constraint Function Values, outdated so does not need to be saved)
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimizer_nlpqlp_data->M; i++)
		fgets(line, sizeof line, fp); 
	fgets(line, sizeof line, fp);  // Space


	// dF (Gradient of the objective function, outdated so does not need to be saved)
	fgets(line, sizeof line, fp);  // Header
	fgets(line, sizeof line, fp); 
	fgets(line, sizeof line, fp);  // Space


	// dG (Gradient of the constraint functions, outdated so does not need to be saved)
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimizer_nlpqlp_data->M; i++)
		fgets(line, sizeof line, fp); 
	fgets(line, sizeof line, fp);  // Space


	// C 
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimizer_nlpqlp_data->NMAX; i++){
		fgets(line, sizeof line, fp);

		strings = line;
		for (j = 0; j < optimizer_nlpqlp_data->NMAX; j++){
			tmpd = strtod(strings, &stringe);
			strings = stringe;
			get_col_Multiarray_d(j,optimizer_nlpqlp_data->C)[i] = tmpd;

		}
	}
	fgets(line, sizeof line, fp);  // Space

	// U 
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimizer_nlpqlp_data->MNN2; i++){
		fgets(line, sizeof line, fp);

		strings = line;
		tmpd = strtod(strings, &stringe);
		strings = stringe;

		optimizer_nlpqlp_data->U->data[i] = tmpd;

	}
	fgets(line, sizeof line, fp);  // Space


	// D
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimizer_nlpqlp_data->NMAX; i++){
		fgets(line, sizeof line, fp);

		strings = line;
		tmpd = strtod(strings, &stringe);
		strings = stringe;

		optimizer_nlpqlp_data->D->data[i] = tmpd;
	}
	fgets(line, sizeof line, fp);  // Space


	// WA
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimizer_nlpqlp_data->LWA; i++){
		fgets(line, sizeof line, fp);

		strings = line;
		tmpd = strtod(strings, &stringe);
		strings = stringe;

		optimizer_nlpqlp_data->WA->data[i] = tmpd;
	}
	fgets(line, sizeof line, fp);  // Space


	// KWA
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimizer_nlpqlp_data->LKWA; i++){
		fgets(line, sizeof line, fp);

		strings = line;
		tmpl = strtol(strings, &stringe, 10);
		strings = stringe;

		optimizer_nlpqlp_data->KWA->data[i] = (int)tmpl;
	}
	fgets(line, sizeof line, fp);  // Space


	// ACTIVE
	fgets(line, sizeof line, fp);  // Header
	for (i = 0; i < optimizer_nlpqlp_data->LACTIV; i++){
		fgets(line, sizeof line, fp);

		if (strstr(line, "T"))
			optimizer_nlpqlp_data->ACTIVE->data[i] = 1;
		else
			optimizer_nlpqlp_data->ACTIVE->data[i] = 0;

	}
	fgets(line, sizeof line, fp);  // Space

	fclose(fp);
}



static void transfer_control_point_data_NLPQLP(struct Optimization_Case *optimization_case, 
	struct Optimizer_NLPQLP_Data* optimizer_nlpqlp_data, char transfer_direction){

	int vec_index = 0;
	int ctrl_pt_index;
	int n_pts;
	struct Multiarray_d* ctrl_pts;
	struct Multiarray_i* ctrl_pts_opt;
	
	if(optimization_case->sim->nurbs_multipatch){
		for(int patch_index=0; patch_index<optimization_case->sim->nurbs_n_patches; patch_index++){
			ctrl_pts = optimization_case->NURBS_Patch_Data[patch_index].control_points;
			ctrl_pts_opt = optimization_case->NURBS_Patch_Data[patch_index].control_points_optimization;
			n_pts = (int)ctrl_pts_opt->extents[0];
			if(transfer_direction == 'o'){
				// Transfer the control point data into the "o"ptimization_case data structure

				for (int i = 0; i < n_pts; i++){
					// Loop over the design control points

					for (int j = 1; j <= 2; j++){
						// Loop over the degrees of freedom for the design point

						if (!get_col_Multiarray_i(j, ctrl_pts_opt)[i])
							continue; // This degree of freedom of the ith control point cannot be moved

						ctrl_pt_index = get_col_Multiarray_i(0, ctrl_pts_opt)[i];
						get_col_Multiarray_d(j-1, ctrl_pts)[ctrl_pt_index] = optimizer_nlpqlp_data->X->data[vec_index];

						vec_index++;
					}
				}


			} else if (transfer_direction == 'n'){
				// Transfer the control point data into the "n"lpqlp data structure

				for (int i = 0; i < n_pts; i++){
					// Loop over the design points

					for (int j = 1; j <= 2; j++){
						// Loop over the degrees of freedom for the design point

						if (!get_col_Multiarray_i(j, ctrl_pts_opt)[i])
							continue; // This degree of freedom of the ith control point cannot be moved

						ctrl_pt_index = get_col_Multiarray_i(0, ctrl_pts_opt)[i];
						optimizer_nlpqlp_data->X->data[vec_index] = get_col_Multiarray_d(j-1, ctrl_pts)[ctrl_pt_index];
					
						vec_index++;
					}
				}
			
			} else {
				EXIT_UNSUPPORTED;
			}
		}
	}
	else{
		ctrl_pts = optimization_case->geo_data.control_points;
		ctrl_pts_opt = optimization_case->geo_data.control_points_optimization;
		n_pts = (int)ctrl_pts_opt->extents[0];
		if(transfer_direction == 'o'){
			// Transfer the control point data into the "o"ptimization_case data structure

			for (int i = 0; i < n_pts; i++){
				// Loop over the design control points

				for (int j = 1; j <= 2; j++){
					// Loop over the degrees of freedom for the design point

					if (!get_col_Multiarray_i(j, ctrl_pts_opt)[i])
						continue; // This degree of freedom of the ith control point cannot be moved

					ctrl_pt_index = get_col_Multiarray_i(0, ctrl_pts_opt)[i];
					get_col_Multiarray_d(j-1, ctrl_pts)[ctrl_pt_index] = optimizer_nlpqlp_data->X->data[vec_index];

					vec_index++;
				}
			}


		} else if (transfer_direction == 'n'){
			// Transfer the control point data into the "n"lpqlp data structure

			for (int i = 0; i < n_pts; i++){
				// Loop over the design points

				for (int j = 1; j <= 2; j++){
					// Loop over the degrees of freedom for the design point

					if (!get_col_Multiarray_i(j, ctrl_pts_opt)[i])
						continue; // This degree of freedom of the ith control point cannot be moved

					ctrl_pt_index = get_col_Multiarray_i(0, ctrl_pts_opt)[i];
					optimizer_nlpqlp_data->X->data[vec_index] = get_col_Multiarray_d(j-1, ctrl_pts)[ctrl_pt_index];
				
					vec_index++;
				}
			}

		} else {
			EXIT_UNSUPPORTED;
		}
	}

}
