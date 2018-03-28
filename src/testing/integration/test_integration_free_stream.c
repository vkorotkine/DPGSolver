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

#include "const_cast.h"
#include "core.h"
#include "math_functions.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"
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
 *  The following conclusions are noted from preliminary testing using a constant exact solution:
 *  - The assumption of a well constructed mesh (\ref CORRECT_NON_CONFORMING_GEOMETRY set to `true`; assumption for
 *    Theorem 2, \cite Kopriva2006) is required for this test to pass.
 *  - The assumption that the geometry order must not be higher than the solution order (Assumption for Theorem 3,
 *    \cite Kopriva2006) is required if the cubature order is limited to 2p. In 2D, using superparametric geometry of
 *    one order higher (p_g = p_s+1) and GL face cubature results in the special case of still obtaining exact
 *    integration (as GL integrates to 2p+1) and consequently satisfying the weak metric identity. Note that using any
 *    higher geometry order requires an increased cubature order for this test to pass.
 *
 *  The conclusions above remained valid on general non-conforming hp-adapted meshes (i.e. non-conforming in both h and
 *  p).
 *
 *  "high-order" free-stream preservation (using an exact solution which is not a constant) is not achieved for the
 *  current implementation, even when the cubature is exact. This implies, for example, that the exact solution given by
 *  a linear profile advected along a constant advection field cannot be **exactly** represented by implementation.
 *  However, it was observed that the infinity norm of the residual for the linear advection equation was converging
 *  asymptotically at a rate of \f$O(h^{p+2})\f$. Noting that the lhs matrix is on the order of \f$O(h^{1})\f$, this
 *  implies that the computed solution update is converging at \f$O(h^{-1}) O(h^{p+2}) = O(h^{p+1})\f$, the optimal
 *  rate, which can be observed when running the convergence order test for the same mesh as used for the free-stream
 *  preservation test.
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

	struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
	const_cast_b(&test_case->copy_initial_rhs,true);

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

	struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
	assert(test_case->copy_initial_rhs == true);

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
