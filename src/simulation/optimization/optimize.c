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

#include "face.h"
#include "volume.h"

#include "macros.h"
#include "multiarray.h"
#include "intrusive.h"
#include "simulation.h"

#include "objective_functions.h"


void optimize(const struct Simulation* sim, const char*const ctrl_file_name){

	/*
	Optimize the the geometry to minimize a specified functional

	Arguments:
		sim = The simulation data structure
		ctrl_file_name = The name of the control file for the  

	Return:
		- 

	*/

	UNUSED(ctrl_file_name);

	printf("Start optimize\n");

	const struct const_Multiarray_d* Cl_Cd = compute_Cl_Cd(sim);

	printf("Cl = %e , Cd = %e \n", Cl_Cd->data[0], Cl_Cd->data[1]);

	//exit(0);

}

