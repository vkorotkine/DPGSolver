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

#include "petscsys.h"
#include "parmetis.h"

#include "simulation.h"
#include "test_base.h"
#include "test_integration.h"
#include "test_integration_convergence_support.h"
//#include "test_support.h"
//#include "test_support_matrix.h"
//#include "test_support_vector.h"
//
//#include "macros.h"
//#include "definitions_core.h"
//#include "definitions_tol.h"
//
//#include "matrix.h"
//#include "vector.h"
//
//#include "math_functions.h"

#include "definitions_adaptation.h" // ADAPT_0
#include "macros.h"
#include "core.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the MPI library.
 *  \return 0 on success.
 *
 *  TO DO:
 *		Figure out which tests we want to perform.
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	assert_condition_message(argc == 3,"Invalid number of input arguments");

	const char*const ctrl_name = argv[1];

	MPI_Init(NULL,NULL);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);



	PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

	struct Integration_Test_Info*const int_test_info = constructor_Integration_Test_Info(ctrl_name); // destructed

	const int p          = int_test_info->p_ref[0],
	          ml         = int_test_info->ml[0],
	          adapt_type = int_test_info->adapt_type;
	destructor_Integration_Test_Info(int_test_info);
	assert(adapt_type == ADAPT_HP);

	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);

	struct Simulation* sim = NULL;

	structor_simulation(&sim,'c',adapt_type,p,ml,0,0,ctrl_name_curr,'r',false); // destructed

	int ParMETIS_Return = METIS_OK;
	//int ParMETIS_Return = ParMETIS_V3_PartMeshKway(elmdist,eptr,eind,elmwgt,wgtflag,numflag,ncon,ncommonnodes,nparts,
	//                                               tpwgts,ubvec,options,edgecut,part,&comm);
	if (ParMETIS_Return != METIS_OK) {
		printf("Parmetis error: %d.\n",ParMETIS_Return);
		EXIT_UNSUPPORTED;
	}

	PetscFinalize();
	MPI_Finalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

