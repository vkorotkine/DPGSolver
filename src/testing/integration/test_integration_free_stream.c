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
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_adaptation.h"
#include "definitions_tol.h"
#include "definitions_visualization.h"

#include "test_base.h"
#include "test_integration.h"

#include "volume_solver.h"

#include "multiarray.h"

#include "core.h"
#include "math_functions.h"
#include "simulation.h"
#include "solve.h"
#include "visualization.h"

// Static function declarations ************************************************************************************* //

/// \brief Check that free-stream preservation is satisfied in each volume.
static void check_free_stream
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for free-stream preservation (\ref test_integration_free_stream.c).
 *  \return 0 on success (when the element-local and global free-stream preservation is satisfied).
 *
 *  Free-stream preservation requires that a constant solution can represented exactly. While this condition is
 *  trivially satisfied on affine meshes, high-order methods have an additional complication of the geometry potentially
 *  introducing flux into the domain when not treated properly.
 *
 *  As a free-stream is represented by a constant solution throughout the domain, the test implemented here checks that
 *  the residual computed from the constant initial solution is equal to zero in each volume of the mesh. For the weak
 *  version of the scheme implemented in the code, this condition is equivalent to the weak satisfaction of the discrete
 *  metric identities presented, for example, by Kopriva (eq. (29), \cite Kopriva2006). Note that the assumption of a
 *  "well-constructed" mesh is also made by Kopriva (which can be tested using \ref test_integration_non_conforming.c).
 *
 *
 *  From preliminary testing, neither the well constructed mesh (\ref CORRECT_NON_CONFORMING_GEOMETRY set to `true`),
 *  nor the use of a geometry order not higher than the solution order (i.e. not superparametric geometry) assumptions
 *  were required for this test to pass; these requirements are discussed by Kopriva (Assumption for Theorem 2 and
 *  Theorem 3, \cite Kopriva2006).
 *
 *  Further, it is noted that computed residual updates were zero (i.e. the test below passed) whenever the exact
 *  solution was of a degree less than or equal to the minimum polynomial degree in all volumes across the mesh. This
 *  can perhaps be interpreted as a high-order free-stream preservation.
 *
 *  \todo update comments after more testing.
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	assert_condition_message(argc == 3,"Invalid number of input arguments");

	const char* petsc_options_name = set_petsc_options_name(argv[2]);
	PetscInitialize(&argc,&argv,petsc_options_name,PETSC_NULL);

	const char* ctrl_name = argv[1];

	struct Test_Info test_info = { .n_warn = 0, };

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);

	const int p  = int_test_info->p_ref[0],
	          ml = int_test_info->ml[0],
	          p_prev  = p-1,
	          ml_prev = ml-1;

	const int adapt_type = int_test_info->adapt_type;
	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);

	struct Simulation* sim = NULL;
	structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr,'r'); // destructed

	solve_for_solution(sim);

	output_visualization(sim,VIS_GEOM_EDGES);
	output_visualization(sim,VIS_SOLUTION);
	check_free_stream(&test_info,sim);

	structor_simulation(&sim,'d',ADAPT_0,p,ml,p_prev,ml_prev,NULL,'r');

	destructor_Integration_Test_Info(int_test_info);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void check_free_stream (struct Test_Info*const test_info, const struct Simulation*const sim)
{
	UNUSED(test_info);
	bool pass = true;
	const double tol = 2e1*EPS;

	double max_norm_rhs = 0.0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;

		const struct Multiarray_d*const rhs = s_vol->rhs;
		const double norm_rhs = norm_Multiarray_d(rhs,"Inf");
		if (norm_rhs > tol) {
			pass = false;
			max_norm_rhs = GSL_MAX(norm_rhs,max_norm_rhs);

			const struct Volume*const vol = (struct Volume*) curr;
			printf("Lack of local free-stream preservation (vol: %d): % .3e \n",vol->index,norm_rhs);
			print_Multiarray_d_tol(rhs,tol);
		}
	}
	if (max_norm_rhs > tol) {
		pass = false;
		printf("Lack of global free-stream preservation: % .3e\n",max_norm_rhs);
	}

	assert_condition(pass);
}
