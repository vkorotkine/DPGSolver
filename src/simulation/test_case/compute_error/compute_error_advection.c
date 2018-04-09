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

#include "compute_error_advection.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "macros.h"
#include "definitions_error.h"

#include "volume.h"
#include "volume_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "compute_error.h"
#include "const_cast.h"
#include "element.h"
#include "intrusive.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Return a statically allocated `char*` holding the specific header for all of the Advection variables.
 *  \return See brief. */
static const char* compute_header_spec_advection_all ( );

/** \brief Return a statically allocated `char*` holding the specific header for all of the Advection variables and the
 *         residual.
 *  \return See brief. */
static const char* compute_header_spec_advection_all_p_rhs ( );

/// \brief Add the computed and exact specified variable to the solution data Multiarrays.
void add_rhs_Error_CE_Data
	(struct Error_CE_Data*const e_ce_d, ///< \ref Error_CE_Data.
	 const struct Simulation*const sim  ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

struct Error_CE* constructor_Error_CE_advection_all (const struct Simulation* sim)
{
	const int n_out = 1;

	struct Error_CE_Helper* e_ce_h = constructor_Error_CE_Helper(sim,n_out);
	e_ce_h->header_spec = compute_header_spec_advection_all();
	const_cast_i(&e_ce_h->error_type,ERROR_STANDARD);

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		e_ce_h->s_vol[0] = (struct Solver_Volume*) curr;
		struct Error_CE_Data* e_ce_d = constructor_Error_CE_Data(e_ce_h,sim); // destructed

		increment_sol_L2(e_ce_h,e_ce_d);
		destructor_Error_CE_Data(e_ce_d);

		update_domain_order(e_ce_h);
	}

	struct Error_CE* error_ce = constructor_Error_CE(e_ce_h,sim); // returned
	destructor_Error_CE_Helper(e_ce_h);

	return error_ce;
}

struct Error_CE* constructor_Error_CE_advection_all_p_rhs (const struct Simulation* sim)
{
	const int n_out = 2;

	struct Error_CE_Helper* e_ce_h = constructor_Error_CE_Helper(sim,n_out);
	e_ce_h->header_spec = compute_header_spec_advection_all_p_rhs();
	const_cast_i(&e_ce_h->error_type,ERROR_STANDARD);

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		e_ce_h->s_vol[0] = (struct Solver_Volume*) curr;
		struct Error_CE_Data* e_ce_d = constructor_Error_CE_Data(e_ce_h,sim); // destructed

		add_rhs_Error_CE_Data(e_ce_d,sim);

		increment_sol_L2(e_ce_h,e_ce_d);
		destructor_Error_CE_Data(e_ce_d);

		update_domain_order(e_ce_h);
	}

	struct Error_CE* error_ce = constructor_Error_CE(e_ce_h,sim); // returned
	destructor_Error_CE_Helper(e_ce_h);

	return error_ce;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const char* compute_header_spec_advection_all ( )
{
	static char header_spec[STRLEN_MAX];
	sprintf(header_spec,"%-14s","$u$");
	return header_spec;
}

static const char* compute_header_spec_advection_all_p_rhs ( )
{
	static char header_spec[STRLEN_MAX];
	sprintf(header_spec,"%-14s%-14s","$u$","$u_{res}$");
	return header_spec;
}

// This function can be used for all pde types if required in future.
void add_rhs_Error_CE_Data (struct Error_CE_Data*const e_ce_d, const struct Simulation*const sim)
{
	struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
	assert(test_case->copy_initial_rhs == true);

	const int n_var = test_case->n_var;

	const ptrdiff_t ext_0 = e_ce_d->sol[0]->extents[0];
	for (int i = 0; i < 2; ++i) {
		struct Multiarray_d*const sol = e_ce_d->sol[i];
//if (i == 0)
//print_Multiarray_d(e_ce_d->rhs[i]);

		const ptrdiff_t ext_1_old = sol->extents[1];
		const ptrdiff_t ext_1_new = ext_1_old+n_var;
		resize_Multiarray_d(sol,sol->order,(ptrdiff_t[]){ext_0,ext_1_new});

		const ptrdiff_t exts_rhs[] = { sol->extents[0], n_var, };
		double*const rhs_ptr = get_col_Multiarray_d(ext_1_old,sol);
		struct Multiarray_d*const sol_rhs =
			constructor_move_Multiarray_d_d('C',sol->order,exts_rhs,true,rhs_ptr); // destructed

		copy_into_Multiarray_d(sol_rhs,(struct const_Multiarray_d*)e_ce_d->rhs[i]);
		sol_rhs->owns_data = false;
		destructor_Multiarray_d(sol_rhs);

		sol->extents[1] = ext_1_new;
	}
}
