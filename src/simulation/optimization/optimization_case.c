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

#include "objective_functions.h"
#include "simulation.h"
#include "volume_solver.h"

#include "test_integration.h"
#include "file_processing.h"

#include "multiarray.h"
#include "complex_multiarray.h"
#include "multiarray_constructors.h"
#include "complex_multiarray_constructors.h"

#include "optimization_case.h"


// Static function declarations ************************************************************************************* //


static void read_geometry_data(struct Optimization_Case* optimization_case);


// Interface functions ********************************************************************************************** //


struct Optimization_Case* constructor_Optimization_Case (const struct Simulation* sim){

	/*
	Create the optimization case structure. Fill it with data relevant to the optimization
	and return the filled structure.

	Arguments:
		sim = The simulation data structure

	Return:
		Optimization_Case data structure with the information for the optimization
	*/

	UNUSED(sim);


	struct Optimization_Case* optimization_case = calloc(1,sizeof *optimization_case); // returned


	// Set function pointers
	optimization_case->objective_function = objective_function_target_Cl;
	optimization_case->objective_function_c = objective_function_target_Cl_c;


	// Read optimization geometry data from .geo file
	read_geometry_data(optimization_case);


	// Create an empty adjoint vector:
	int RHS_size_ex_0;
	struct Solver_Volume* s_vol = (struct Solver_Volume*)sim->volumes->first;
	int P = s_vol->p_ref,  // P value for all elements
		nv = (int)sim->n_v,  // number of volumes  
		neq = NEQ_EULER,  // number of equations
		d = DIM;  // Dimension
	RHS_size_ex_0 = nv * (int)(pow(P+1, d)) * neq;
	optimization_case->Chi = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){RHS_size_ex_0,1});


	return optimization_case;

}


void destructor_Optimization_Case (struct Optimization_Case* optimization_case){

	/*
	Free the optimization_case structure

	Arguments:
		optimization_case = The optimization_case structure

	Return:
		-
	*/

	// Destroy the geometry data
	// TODO: Implement this portion

	// Destroy the simulation. Note that the p and ml values are not used when 
	// the solution is being destructed
	struct Simulation* sim_c = optimization_case->sim_c;
	structor_simulation(&sim_c,'d',ADAPT_0,0,0,0,0,NULL,'c',false);

	destructor_Multiarray_d(optimization_case->Chi);

	free((void*)optimization_case);

}


// Static functions ************************************************************************************************* //

static void read_geometry_data(struct Optimization_Case* optimization_case){

	/*
	Read the data from the geometry_parameters.geo file. For now, only a NURBS patch
	case can be read for the optimization. The .geo file should contain a section
	with the optimization control points.

	Arguments:
		optimization_case = The data structure holding the optimization data

	Return:
		-
	*/

	// For now, 2D patches should be able to be read successfully so
	// we should be able to find
	//	- P (order in the xi direction)
	// 	- Q (order in the eta direction)
	//	- knots_xi
	//  - knots_eta
	// 	- Control Points and Weights
	//  - Control Point Connectivity
	// 	- Optimization_Point_Connectivity

	const int count_to_find = 7;
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
			struct Multiarray_d* control_points_and_weights = 
				constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_control_pts,DIM+1}); // saved

			double *x = get_col_Multiarray_d(0, control_points_and_weights),
				   *y = get_col_Multiarray_d(1, control_points_and_weights),
				   *w = get_col_Multiarray_d(2, control_points_and_weights);

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
			optimization_case->geo_data.control_points_and_weights   = control_points_and_weights;
			optimization_case->geo_data.control_points_and_weights_c = 
				constructor_copy_Multiarray_c_Multiarray_d(control_points_and_weights);
		}

		// Read the Control Point Connecitivity Data
		if (strstr(line,"Control_Point_Connectivity")){

			// Get the number of xi and eta points
			count_found++;  // Connectivity information was found
			int extents[2] = {0,0};
			read_skip_const_i_1(line,1,extents,2);

			// Multiarray structure to hold the connectivity data
			struct Multiarray_i* control_point_connectivity = 
				constructor_empty_Multiarray_i('R',2,(ptrdiff_t[]){extents[0], extents[1]});

			// Read in the connectivity data line by line
			for (i = 0; i < extents[0]; i++){
				fgets(line,sizeof(line),input_file);
				read_skip_i_1(line,0, get_row_Multiarray_i(i, control_point_connectivity), extents[1]);
			}

			optimization_case->geo_data.control_point_connectivity = control_point_connectivity;

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
	}

	fclose(input_file);

	assert(count_found == count_to_find);

	/*
		// Test that the information is read properly
	printf("P = %d, Q = %d \n", optimization_case->geo_data.P, optimization_case->geo_data.Q);
	print_Multiarray_d(optimization_case->geo_data.knots_xi);
	print_Multiarray_d(optimization_case->geo_data.knots_eta);		
	print_Multiarray_d(optimization_case->geo_data.control_points_and_weights);
	print_Multiarray_c(optimization_case->geo_data.control_points_and_weights_c);
	print_Multiarray_i(optimization_case->geo_data.control_point_connectivity);
	print_Multiarray_i(optimization_case->geo_data.control_points_optimization);
	
	exit(0);
	*/

}





