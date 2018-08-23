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
#include "petscsys.h"

#include "macros.h"

#include "test_base.h"
#include "test_integration.h"
#include "test_integration_convergence_support.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for restarted solutions (\ref test_integration_restart.c).
 *  \return 0 on success.
 *  \note The appropriate restart file must be created using \ref test_integration_restart_output.c before this test can
 *        run (as the restart file will otherwise not be present).
 *
 *  Requires a test case with a non-constant specified exact solution. This test loops over the all degrees and mesh
 *  levels and computes the error between the restarted and exact solution.
 *
 *  \note It was found that using a restart file which was more refined sometimes resulted in larger errors and
 *        stagnation of errors well above the machine precision. This occured regardless of: the restarted solution
 *        being interpolated or L2-projected from the restart file solution, the approximate nearest neighbor algorithm
 *        using `int` or `long long int` for the integer coordinate representation. It was also noted that using a
 *        restart file with triangular elements led to significantly higher error than quad elements for the Euler
 *        supersonic_vortex case; the problem thus may potentially be coming from the parametric domain mapping.
 *  \todo Investigate why this is occuring.
 *
 *  \warning The test is failing for p > 2 when using triangular elements. Only entropy was converging optimally in a
 *  previously tested euler test case while all other variables were converging at sub-optimal rates. The problem was
 *  not observed when using only quad meshes (i.e. convergence was optimal). Is this related to the loss of 2p cubature
 *  strength for triangular elements start at p = 3? Requires investigation.
 *
 *  The test passes if the convergence is optimal (O(h^{p+1}).
 *
 *  Please see comments in \ref restart_T.c regarding potential superconvergence of functionals depending on the
 *  projection method for the solution from the restart file.
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
	assert_condition_message(argc == 2,"Invalid number of input arguments");

	run_convergence_order_study(argc,argv,CONV_STUDY_RESTART);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
