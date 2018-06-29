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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"

#include "definitions_core.h"
#include "definitions_adaptation.h"
#include "definitions_core.h"

//#include "objective_functions.h"
#include "simulation.h"
#include "volume_solver.h"
#include "face_solver.h"
#include "intrusive.h"
#include "definitions_intrusive.h"
#include "face.h"
#include "volume.h"

#include "test_integration.h"
#include "file_processing.h"

#include "multiarray.h"
#include "multiarray_constructors.h"

#include "optimization_case.h"


// Static function declarations ************************************************************************************* //

/** \brief 	Read the data from the geometry_parameters.geo file. For now, only a NURBS patch
 * case can be read for the optimization. The .geo file should contain a section
 * with the optimization control points and their upper and lower limits. 
 * 
 * NOTE: This function is very similar to read_NURBS_data in geometry_parametric however it
 * has additional data to read regarding the optimization information. Look to perhaps 
 * combine these functions.
 */
static void read_geometry_data(
	struct Optimization_Case* optimization_case ///< The data structure holding the optimization data
	);


// Interface functions ********************************************************************************************** //


struct Optimization_Case* constructor_Optimization_Case (struct Simulation* sim){

	struct Optimization_Case* optimization_case = calloc(1,sizeof *optimization_case); // returned


	// Allocate the memory in the Optimization_Case structure except for those used by the optimizers
	// and adjoint solve (Adjoint Data Structures, Sensitivity Data Structures, Optimization Minimizer Structures
	// NLPQLP Data Structures). Load the information for the optimization case into the relevant design variables


	// 1) Geometry Data Structures
	read_geometry_data(optimization_case);


	// 2) num_design_pts_dofs
	int num_design_pt_dofs = 0;

	struct Multiarray_i* ctrl_pts_opt = optimization_case->geo_data.control_points_optimization;
	int n_pts = (int)ctrl_pts_opt->extents[0];
	
	for (int i = 0; i < n_pts; i++){
		// Loop over the design points

		for (int j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point

			if (get_col_Multiarray_i(j, ctrl_pts_opt)[i])
				num_design_pt_dofs++;
		}
	}	

	optimization_case->num_design_pts_dofs = num_design_pt_dofs;


	// 3) objective_function and objective_function_c
	// Set function pointers
	// MSB: TODO: Set these based on optimization.data file
	optimization_case->objective_function = objective_function_target_Cl;
	optimization_case->objective_function_c = objective_function_target_Cl_c;


	// 4) Sim objects (sim, sim_c)

	optimization_case->sim = sim;

	// Create the complex simulation object

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(sim->ctrl_name); // free
	
	const int* p_ref  = int_test_info->p_ref,
	         * ml_ref = int_test_info->ml;

	int ml_prev = ml_ref[0]-1,
	    ml 		= ml_ref[0],
    	p_prev  = p_ref[0]-1,
    	p 		= p_ref[0];

	const char type_rc = 'c';
	bool ignore_static = false;
	const int adapt_type = int_test_info->adapt_type;

	struct Simulation* sim_c = NULL;
	structor_simulation(&sim_c, 'c', adapt_type, p, ml, p_prev, ml_prev, sim->ctrl_name, type_rc, ignore_static);
	
	optimization_case->sim_c = sim_c;


	// Copy the data from the real to the complex sim object
	// \todo MSB: For the copy members r_to_c_Solver_Face function, it seems that the nf_fc 
	// data structures needs to be set here to operate correctly
	struct Intrusive_Link* curr   = sim->faces->first;
	struct Intrusive_Link* curr_c = sim_c->faces->first;

	while(true){

		struct Solver_Face* s_face 		= (struct Solver_Face*) curr;
		struct Solver_Face_c* s_face_c 	= (struct Solver_Face_c*) curr_c;
		
		// TODO:
		// Needed in order to do the r_to_c conversion. Find a way to fix this issue
		s_face->nf_fc = (const struct const_Multiarray_d*)constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1,1});
		s_face_c->nf_fc = (const struct const_Multiarray_c*)constructor_empty_Multiarray_c('C',2,(ptrdiff_t[]){1,1});

		curr 	= curr->next;
		curr_c 	= curr_c->next;

		if(!curr || !curr_c)
			break;
	}

	copy_data_r_to_c_sim(optimization_case->sim, optimization_case->sim_c);


	// 5) Optimizer
	// MSB: TODO: Read from optimization.data file
	//strcpy(optimization_case->optimizer_spec, "LINE_SEARCH_BFGS");
	strcpy(optimization_case->optimizer_spec, "NLPQLP");

	// 6) Optimization Type
	// MSB: TODO: Read from optimization.data file
	strcpy(optimization_case->optimization_type_spec, "TARGET_CL");

	
	destructor_Integration_Test_Info(int_test_info);

	return optimization_case;

}


void destructor_Optimization_Case (struct Optimization_Case* optimization_case){

	// 1) Geometry Data Structures
	destructor_Multiarray_d(optimization_case->geo_data.knots_xi);
	destructor_Multiarray_d(optimization_case->geo_data.knots_eta);
	destructor_Multiarray_d(optimization_case->geo_data.control_points);
	destructor_Multiarray_c(optimization_case->geo_data.control_points_c);
	destructor_Multiarray_d(optimization_case->geo_data.control_weights);
	destructor_Multiarray_i(optimization_case->geo_data.control_pt_wt_connectivity);
	destructor_Multiarray_i(optimization_case->geo_data.control_points_optimization);
	destructor_Multiarray_d(optimization_case->geo_data.control_points_optimization_lims);


	// 2) Sim_c object
	// Destroy the simulation. Note that the p and ml values are not used when 
	// the solution is being destructed

	// NOTE: Real sim object not cleared
	struct Simulation* sim_c = optimization_case->sim_c;
	structor_simulation(&sim_c,'d',ADAPT_0,0,0,0,0,NULL,'c',false);

	free((void*)optimization_case);

}


void copy_data_r_to_c_sim(struct Simulation* sim, struct Simulation* sim_c){

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
		
		copy_members_r_to_c_Solver_Face((struct Solver_Face_c*const)s_face_c, 
			(const struct Solver_Face*const)s_face, sim);

		curr 	= curr->next;
		curr_c 	= curr_c->next;

		if(!curr || !curr_c)
			break;
	}
}


// Static functions ************************************************************************************************* //

static void read_geometry_data(struct Optimization_Case* optimization_case){

	// For now, 2D patches should be able to be read successfully so
	// we should be able to find
	//	- P (order in the xi direction)
	// 	- Q (order in the eta direction)
	//	- knots_xi
	//  - knots_eta
	// 	- Control Points and Weights
	//  - Control Point Connectivity
	// 	- Optimization_Point_Connectivity
	// 	- Optimization_Point_Limit

	const int count_to_find = 8;
	assert(DIM == 2);

	// Get the file pointer to the geometry file
	FILE* input_file = fopen_input('g',NULL,NULL); // closed

	int count_found = 0;
	char line[STRLEN_MAX];

	int i, len_knots_xi, len_knots_eta;

	// Read in the information from the file
	while (fgets(line,sizeof(line),input_file)) {

		// Get the P and Q information
		read_skip_string_count_i("P(xi_order)", &count_found, line, &optimization_case->geo_data.P);
		read_skip_string_count_i("Q(eta_order)", &count_found, line, &optimization_case->geo_data.Q);

		// The knot vectors in the eta and xi directions:
		if (strstr(line, "knots_xi")){
			read_skip_string_count_i("knots_xi", &count_found, line, &len_knots_xi);
				
			struct Multiarray_d* knots_xi = 
					constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){len_knots_xi,1});
			
			for (i = 0; i < len_knots_xi; i++){
				fgets(line,sizeof(line),input_file);  // read the next line
				read_skip_d_1(line, 0, &get_col_Multiarray_d(0, knots_xi)[i], 1);
			}

			optimization_case->geo_data.knots_xi = knots_xi;
		}

		if (strstr(line, "knots_eta")){
			read_skip_string_count_i("knots_eta", &count_found, line, &len_knots_eta);
			
			struct Multiarray_d* knots_eta = 
					constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){len_knots_eta,1});
			
			for (i = 0; i < len_knots_eta; i++){
				fgets(line,sizeof(line),input_file);  // read the next line
				read_skip_d_1(line, 0, &get_col_Multiarray_d(0, knots_eta)[i], 1);
			}
			
			optimization_case->geo_data.knots_eta = knots_eta;

		}

		// Read the Control Point Data
		if (strstr(line,"Control_Point_Data")){
			
			// Get the number of pts
			int num_control_pts;
			read_skip_string_count_i("Control_Point_Data", &count_found, line, &num_control_pts);

			// Create the multiarray structure to hold the control point data
			struct Multiarray_d *control_points  = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_control_pts,DIM}); // saved
			struct Multiarray_d *control_weights = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_control_pts,1}); // saved

			double *x = get_col_Multiarray_d(0, control_points),
				   *y = get_col_Multiarray_d(1, control_points),
				   *w = get_col_Multiarray_d(0, control_weights);

			// Read the pt data line by line
			for (i = 0; i < num_control_pts; i++){
				fgets(line,sizeof(line),input_file);
				
				double pt[3] = {0.0, 0.0, 0.0};
				read_skip_d_1(line,0, pt, DIM+1);

				x[i] = pt[0];
				y[i] = pt[1];
				w[i] = pt[2];
			}

			// Save the control pt data and store a complex copy
			optimization_case->geo_data.control_points   = control_points;
			optimization_case->geo_data.control_points_c = 
				constructor_copy_Multiarray_c_Multiarray_d(control_points);

			optimization_case->geo_data.control_weights  = control_weights;
		}

		// Read the Control Point Connecitivity Data
		if (strstr(line,"Control_Point_Connectivity")){

			// Get the number of xi and eta points
			count_found++;  // Connectivity information was found
			int extents[2] = {0,0};
			read_skip_const_i_1(line,1,extents,2);

			// Multiarray structure to hold the connectivity data
			struct Multiarray_i* control_pt_wt_connectivity = 
				constructor_empty_Multiarray_i('R',2,(ptrdiff_t[]){extents[0], extents[1]});

			// Read in the connectivity data line by line
			for (i = 0; i < extents[0]; i++){
				fgets(line,sizeof(line),input_file);
				read_skip_i_1(line,0, get_row_Multiarray_i(i, control_pt_wt_connectivity), extents[1]);
			}

			optimization_case->geo_data.control_pt_wt_connectivity = control_pt_wt_connectivity;

		}

		// Read the Optimization Point Data
		if (strstr(line,"Optimization_Point_Connectivity")){
			
			int num_opt_pts;
			read_skip_string_count_i("Optimization_Point_Connectivity", &count_found, line, &num_opt_pts);

			// Multiarray structure to hold the optimization control point data.
			struct Multiarray_i* optimization_pts = 
				constructor_empty_Multiarray_i('C',2,(ptrdiff_t[]){num_opt_pts,DIM+1}); // saved

			int *dof_i 		= get_col_Multiarray_i(0, optimization_pts),  // connectivity integer
				*x_dof_i 	= get_col_Multiarray_i(1, optimization_pts),  // x direction dof bool
				*y_dof_i 	= get_col_Multiarray_i(2, optimization_pts);  // y direction dof bool

			// Read the pt data line by line
			for (i = 0; i < num_opt_pts; i++){
				fgets(line,sizeof(line),input_file);
				
				int line_data[3] = {0, 0, 0};
				read_skip_i_1(line, 0, line_data, DIM+1);

				dof_i  [i] = line_data[0];
				x_dof_i[i] = line_data[1];
				y_dof_i[i] = line_data[2];
			}

			optimization_case->geo_data.control_points_optimization = optimization_pts;
		}

		// Read the Optimization Point Limits Data
		if (strstr(line,"Optimization_Point_Limit")){
			
			int num_opt_pts;
			read_skip_string_count_i("Optimization_Point_Limit", &count_found, line, &num_opt_pts);

			// Multiarray structure to hold the optimization control point data.
			struct Multiarray_d* optimization_pts_lims = 
				constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_opt_pts,2}); // saved

			double 	*low_lim 	= get_col_Multiarray_d(0, optimization_pts_lims),
					*high_lim 	= get_col_Multiarray_d(1, optimization_pts_lims);

			// Read the pt data line by line
			for (i = 0; i < num_opt_pts; i++){
				fgets(line,sizeof(line),input_file);
				
				double line_data[3] = {0, 0, 0};
				read_skip_d_1(line, 0, line_data, 3);

				low_lim[i]  = line_data[1];
				high_lim[i] = line_data[2];
			}

			optimization_case->geo_data.control_points_optimization_lims = optimization_pts_lims;
		}
	}

	fclose(input_file);

	assert(count_found == count_to_find);

	// // Test that the information is read properly
	// printf("P = %d, Q = %d \n", optimization_case->geo_data.P, optimization_case->geo_data.Q);
	// print_Multiarray_d(optimization_case->geo_data.knots_xi);
	// print_Multiarray_d(optimization_case->geo_data.knots_eta);		
	// print_Multiarray_d(optimization_case->geo_data.control_points);
	// print_Multiarray_c(optimization_case->geo_data.control_points_c);
	// print_Multiarray_d(optimization_case->geo_data.control_weights);
	// print_Multiarray_i(optimization_case->geo_data.control_pt_wt_connectivity);
	// print_Multiarray_i(optimization_case->geo_data.control_points_optimization);
	// print_Multiarray_d(optimization_case->geo_data.control_points_optimization_lims);
	// exit(0);
}





