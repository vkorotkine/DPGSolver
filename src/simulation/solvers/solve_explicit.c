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

#include "solve_explicit.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_test_case.h"
#include "definitions_tol.h"

#include "computational_elements.h"
#include "volume_solver.h"
#include "volume_solver_dg.h"

#include "multiarray.h"

#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Function pointer to a time-stepping function.
 *  \return The absolute value of the maximum rhs at the current time.
 *  \param dt  The time step.
 *  \param sim \ref Simulation.
 */
typedef double (*time_step_fptr)
	(const double dt,
	 const struct Simulation* sim
	);

/** \brief Set the function pointer to the appropriate time-stepping function.
 *  \return See brief. */
static time_step_fptr set_time_step
	(const struct Simulation* sim ///< \ref Simulation
	);

/// \brief Display the solver progress.
static void display_progress
	(const struct Test_Case* test_case, ///< \ref Test_Case_T.
	 const int t_step,                  ///< The current time step.
	 const double max_rhs,              ///< The current maximum value of the rhs term.
	 const double max_rhs0              ///< The initial maximum value of the rhs term.
	);

/** \brief Check the exit conditions.
 *  \return `true` if exit conditions are satisfied, `false` otherwise. */
static bool check_exit
	(const struct Test_Case* test_case, ///< \ref Test_Case_T.
	 const double max_rhs,              ///< The current maximum value of the rhs term.
	 const double max_rhs0              ///< The initial maximum value of the rhs term.
	);

// Interface functions ********************************************************************************************** //

void solve_explicit (struct Simulation* sim)
{
	assert(sim->method == METHOD_DG); // Can be made flexible in future.

	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	test_case->solver_method_curr = 'e';

	test_case->time = 0.0;
	set_initial_solution(sim);

	constructor_derived_Elements(sim,IL_ELEMENT_SOLVER_DG);       // destructed
	constructor_derived_computational_elements(sim,IL_SOLVER_DG); // destructed

	time_step_fptr time_step = set_time_step(sim);

	const double time_final = test_case->time_final;
	double dt = test_case->dt;
	assert(time_final >= 0.0);

	double max_rhs0 = 0.0;
	for (int t_step = 0; ; ++t_step) {
		if (test_case->time + dt > time_final-1e3*EPS)
			dt = time_final - test_case->time;
		test_case->time += dt;

		const double max_rhs = time_step(dt,sim);
		if (t_step == 0) {
			max_rhs0 = max_rhs;
			if (test_case->copy_initial_rhs)
				copy_rhs(sim,NULL);
		}

		display_progress(test_case,t_step,max_rhs,max_rhs0);
		if (check_exit(test_case,max_rhs,max_rhs0))
			break;
	}

	destructor_derived_computational_elements(sim,IL_SOLVER);
	destructor_derived_Elements(sim,IL_ELEMENT_SOLVER);
	test_case->solver_method_curr = 0;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Version of \ref time_step_fptr using the forward Euler method.
 *  \return See \ref time_step_fptr. */
static double time_step_euler
	(const double dt,             ///< The time step.
	 const struct Simulation* sim ///< Defined for \ref time_step_fptr.
	);

/** \brief Version of \ref time_step_fptr using the strong stability preserving Runge-Kutta 3-stage, 3rd order method.
 *  \return See \ref time_step_fptr.
 *
 *  Reference: eq. (4.2) \cite Gottlieb2001.
 */
static double time_step_ssp_rk_33
	(const double dt,             ///< The time step.
	 const struct Simulation* sim ///< Defined for \ref time_step_fptr.
	);

/** \brief Version of \ref time_step_fptr using the low-storage Runge-Kutta 5-stage, 4rd order method.
 *  \return See \ref time_step_fptr.
 *
 *  Reference: Solution 3, Table 1, rational form \cite Carpenter1994.
 */
static double time_step_ls_rk_54
	(const double dt,             ///< The time step.
	 const struct Simulation* sim ///< Defined for \ref time_step_fptr.
	);

static time_step_fptr set_time_step (const struct Simulation* sim)
{
	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	assert((test_case->solver_proc == SOLVER_E) || (test_case->solver_proc == SOLVER_EI));

	switch (test_case->solver_type_e) {
	case SOLVER_E_EULER:
		return time_step_euler;
		break;
	case SOLVER_E_SSP_RK_33:
		return time_step_ssp_rk_33;
		break;
	case SOLVER_E_LS_RK_54:
		return time_step_ls_rk_54;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->solver_type_e);
		break;
	}
}

static void display_progress
	(const struct Test_Case* test_case, const int t_step, const double max_rhs, const double max_rhs0)
{
	if (!test_case->display_progress)
		return;

	switch (test_case->solver_proc) {
	case SOLVER_E:
		printf("Complete: % 7.2f%%, tstep: %8d, maxRHS: % .3e\n",
		       100*(test_case->time)/(test_case->time_final),t_step,max_rhs);
		break;
	case SOLVER_EI:
		printf("Exit Conditions (tol, ratio): % .3e (% .3e), % .3e (% .3e), tstep: %8d\n",
		       max_rhs,test_case->exit_tol_e,max_rhs/max_rhs0,test_case->exit_ratio_e,t_step);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->solver_proc);
		break;
	}
}

static bool check_exit (const struct Test_Case* test_case, const double max_rhs, const double max_rhs0)
{
	bool exit_now = false;

	if (max_rhs < test_case->exit_tol_e) {
		printf("Complete: max_rhs is below the exit tolerance.\n");
		exit_now = true;
	}

	switch (test_case->solver_proc) {
	case SOLVER_E:
		if (test_case->time == test_case->time_final)
			exit_now = true;
		break;
	case SOLVER_EI:
		if (max_rhs < test_case->exit_tol_e) {
			printf("Complete: max_rhs is below the exit tolerance.\n");
			exit_now = true;
		}

		if (max_rhs/max_rhs0 < test_case->exit_ratio_e) {
			printf("Complete: max_rhs dropped by % .2e orders.\n",log10(max_rhs0/max_rhs));
			exit_now = true;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->solver_proc);
		break;
	}
	return exit_now;
}

// Level 1 ********************************************************************************************************** //

static double time_step_euler (const double dt, const struct Simulation* sim)
{
	assert(sim->method == METHOD_DG);
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DG);
	assert(sim->faces->name   == IL_FACE_SOLVER_DG);

	const double max_rhs = compute_rhs(sim);
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		// s_vol holds the data for the given volume
		// s_vol_dg holds the dg related data, such as the RHS for the given volume
		struct Solver_Volume*    s_vol    = (struct Solver_Volume*) curr;
		struct DG_Solver_Volume* s_vol_dg = (struct DG_Solver_Volume*) curr;

		struct Multiarray_d* sol_coef = s_vol->sol_coef,
		                   * rhs      = s_vol_dg->rhs;

		double* data_s   = sol_coef->data,
		      * data_rhs = rhs->data;

		const ptrdiff_t i_max = compute_size(sol_coef->order,sol_coef->extents);
		assert(i_max == compute_size(rhs->order,rhs->extents));

		for (ptrdiff_t i = 0; i < i_max; ++i)
			data_s[i] += dt*data_rhs[i];
		enforce_positivity_highorder(s_vol,sim);
	}
	return max_rhs;
}

static double time_step_ssp_rk_33 (const double dt, const struct Simulation* sim)
{
	assert(sim->method == METHOD_DG);
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DG);
	assert(sim->faces->name   == IL_FACE_SOLVER_DG);

	double max_rhs = 0.0;
	for (int rk = 0; rk < 3; rk++) {
		max_rhs = compute_rhs(sim);
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Solver_Volume*    s_vol    = (struct Solver_Volume*) curr;
			struct DG_Solver_Volume* s_vol_dg = (struct DG_Solver_Volume*) curr;

			struct Multiarray_d* sol_coef   = s_vol->sol_coef,
			                   * sol_coef_p = s_vol_dg->sol_coef_p,
			                   * rhs        = s_vol_dg->rhs;

			double* data_s   = sol_coef->data,
			      * data_sp  = sol_coef_p->data,
			      * data_rhs = rhs->data;

			const ptrdiff_t i_max = compute_size(sol_coef->order,sol_coef->extents);
			assert(i_max == compute_size(sol_coef_p->order,sol_coef_p->extents));
			assert(i_max == compute_size(rhs->order,rhs->extents));

			switch (rk) {
			case 0:
				for (ptrdiff_t i = 0; i < i_max; ++i) {
					data_sp[i]  = data_s[i];
					data_s[i]  += dt*data_rhs[i];
				}
				break;
			case 1:
				for (ptrdiff_t i = 0; i < i_max; ++i)
					data_s[i] = (1.0/4.0)*(3.0*data_sp[i] + data_s[i] + dt*data_rhs[i]);
				break;
			case 2:
				for (ptrdiff_t i = 0; i < i_max; ++i)
					data_s[i] = (1.0/3.0)*(data_sp[i] + 2.0*data_s[i] + 2.0*dt*data_rhs[i]);
				break;
			default:
				EXIT_ERROR("Unsupported: %d\n",rk);
				break;
			}
			enforce_positivity_highorder(s_vol,sim);
		}
	}

	return max_rhs;
}

static double time_step_ls_rk_54 (const double dt, const struct Simulation* sim)
{
	assert(sim->method == METHOD_DG);
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DG);
	assert(sim->faces->name   == IL_FACE_SOLVER_DG);

	static const double rk4a[] =
		{  0.0,                              -5673018057730.0/13575370590870.0, -2404267990393.0/20167466952380.0,
		  -3550918686646.0/20915011793850.0, -1275806237668.0/84257045769900.0, };
	static const double rk4b[] =
		{  1432997174477.0/95750804417550.0,  5161836677717.0/13612068292357.0, 17201463215490.0/20902069494980.0,
		   3134564353537.0/44814673103380.0,  2277821191437.0/14882151754819.0, };
	static const double rk4c[] =
		{  0.0,                               1432997174477.0/95750804417550.0,  2526269341429.0/68203639628960.0,
		   2006345519317.0/32243100637760.0,  2802321613138.0/29243179262510.0, };
	UNUSED(rk4c);

	double max_rhs = 0.0;
	for (int rk = 0; rk < 5; rk++) {
		max_rhs = compute_rhs(sim);
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Solver_Volume*    s_vol    = (struct Solver_Volume*) curr;
			struct DG_Solver_Volume* s_vol_dg = (struct DG_Solver_Volume*) curr;

			struct Multiarray_d* sol_coef   = s_vol->sol_coef,
			                   * sol_coef_p = s_vol_dg->sol_coef_p,
			                   * rhs        = s_vol_dg->rhs;

			double* data_s   = sol_coef->data,
			      * data_sp  = sol_coef_p->data,
			      * data_rhs = rhs->data;

			const ptrdiff_t i_max = compute_size(sol_coef->order,sol_coef->extents);
			assert(i_max == compute_size(sol_coef_p->order,sol_coef_p->extents));
			assert(i_max == compute_size(rhs->order,rhs->extents));

			for (ptrdiff_t i = 0; i < i_max; ++i) {
				data_sp[i] *= rk4a[rk];
				data_sp[i] += dt*data_rhs[i];
				data_s[i]  += rk4b[rk]*data_sp[i];
			}
			enforce_positivity_highorder(s_vol,sim);
		}
	}

	return max_rhs;
}
