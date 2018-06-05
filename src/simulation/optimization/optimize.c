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

#include "macros.h"

#include "multiarray.h"
#include "complex_multiarray.h"
#include "multiarray_constructors.h"
#include "complex_multiarray_constructors.h"

#include "intrusive.h"
#include "definitions_intrusive.h"

#include "simulation.h"

#include "volume_solver.h"
#include "face_solver.h"
#include "test_complex_face_solver.h"
#include "test_complex_volume_solver.h"

#include "optimization_case.h"

#include "adjoint.h"
#include "sensitivities.h"


// Static function declarations ************************************************************************************* //

static struct Optimization_Case* setup_optimization(struct Simulation* sim);

static void copy_data_r_to_c_sim(struct Simulation* sim, struct Simulation* sim_c);


// Interface functions ********************************************************************************************** //

void optimize(struct Simulation* sim){

	/*
	Optimize the the geometry to minimize a specified functional

	Arguments:
		sim = The simulation data structure

	Return:
		- 

	*/

	// Preprocessing:
	struct Optimization_Case *optimization_case = setup_optimization(sim);
	copy_data_r_to_c_sim(optimization_case->sim, optimization_case->sim_c);


	// Optimization Routine

	setup_adjoint(optimization_case);
	solve_adjoint(optimization_case);

	compute_sensitivities(optimization_case);

	printf("dI_dXp: \n");
	print_Multiarray_d(optimization_case->dI_dXp);


	double I_val = optimization_case->objective_function(optimization_case->sim);
	double complex I_val_c = optimization_case->objective_function_c(optimization_case->sim_c);
	
	printf("I   : %e \n", I_val);
	printf("I_c : %e + %e * i\n", creal(I_val_c), cimag(I_val_c));




	// Clear allocated structures:
	destructor_Optimization_Case(optimization_case);

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








