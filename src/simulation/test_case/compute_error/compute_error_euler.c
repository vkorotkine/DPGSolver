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

#include "compute_error_euler.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "macros.h"
#include "definitions_core.h"

#include "volume.h"
#include "volume_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "compute_error.h"
#include "const_cast.h"
#include "element.h"
#include "intrusive.h"
#include "simulation.h"
#include "solution_euler.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Return a statically allocated `char*` holding the specific header for all of the Euler variables.
 *  \return See brief. */
static const char* compute_header_spec_euler_all
	();

/** \brief Return a statically allocated `char*` holding the specific header for the entropy variable.
 *  \return See brief. */
static const char* compute_header_spec_euler_entropy
	();

// Interface functions ********************************************************************************************** //

struct Error_CE* constructor_Error_CE_euler_all (const struct Simulation* sim)
{
	const int n_out = DIM+2+1;

	struct Error_CE_Helper* e_ce_h = constructor_Error_CE_Helper(sim,n_out);
	e_ce_h->header_spec = compute_header_spec_euler_all();

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		e_ce_h->s_vol[0] = (struct Solver_Volume*) curr;
		struct Error_CE_Data* e_ce_d = constructor_Error_CE_Data(e_ce_h,sim); // destructed

		for (int i = 0; i < 2; ++i)
			convert_variables(e_ce_d->sol[i],'c','p');
		add_euler_variable_Error_CE_Data('s',e_ce_d,sim);

		increment_sol_L2(e_ce_h,e_ce_d);
		destructor_Error_CE_Data(e_ce_d);

		update_domain_order(e_ce_h);
	}

	struct Error_CE* error_ce = constructor_Error_CE(e_ce_h,sim); // returned
	destructor_Error_CE_Helper(e_ce_h);

	return error_ce;
}

struct Error_CE* constructor_Error_CE_euler_entropy (const struct Simulation* sim)
{
	const int n_out = 1;

	struct Error_CE_Helper* e_ce_h = constructor_Error_CE_Helper(sim,n_out);
	e_ce_h->header_spec = compute_header_spec_euler_entropy();

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		e_ce_h->s_vol[0] = (struct Solver_Volume*) curr;
		struct Error_CE_Data* e_ce_d = constructor_Error_CE_Data(e_ce_h,sim); // destructed

		for (int i = 0; i < 2; ++i)
			convert_variables(e_ce_d->sol[i],'c','p');
		add_euler_variable_Error_CE_Data('s',e_ce_d,sim);
		for (int i = 0; i < 2; ++i) {
			for (int vr = 0; vr < DIM+1; ++vr) // Remove all except entropy.
				remove_col_Multiarray_d(0,e_ce_d->sol[i]);
		}

		increment_sol_L2(e_ce_h,e_ce_d);
		destructor_Error_CE_Data(e_ce_d);

		update_domain_order(e_ce_h);
	}

	struct Error_CE* error_ce = constructor_Error_CE(e_ce_h,sim); // returned
	destructor_Error_CE_Helper(e_ce_h);

	return error_ce;
}

void add_euler_variable_Error_CE_Data
	(const char var_type, struct Error_CE_Data*const e_ce_d, const struct Simulation*const sim)
{
	UNUSED(sim);
	const ptrdiff_t ext_0 = e_ce_d->sol[0]->extents[0];
	struct Multiarray_d*const var = constructor_move_Multiarray_d_d('C',2,(ptrdiff_t[]){ext_0,1},false,NULL); // dest.
	for (int i = 0; i < 2; ++i) {
		struct Multiarray_d* sol = e_ce_d->sol[i];

		const ptrdiff_t ext_1_new = sol->extents[1]+1;
		resize_Multiarray_d(sol,sol->order,(ptrdiff_t[]){ext_0,ext_1_new});

		sol->extents[1] = DIM+2;
		var->data = get_col_Multiarray_d(ext_1_new-1,sol);
		switch (var_type) {
			case 's': compute_entropy(var,(const struct const_Multiarray_d*)sol,'p'); break;
			case 't': compute_temperature(var,(const struct const_Multiarray_d*)sol,'p'); break;
			default:  EXIT_ERROR("Unsupported: %c.",var_type); break;
		}
		sol->extents[1] = ext_1_new;
	}
	destructor_Multiarray_d(var);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const char* compute_header_spec_euler_all ( )
{
	static char header_spec[STRLEN_MAX];

	int index = sprintf(header_spec,"%-14s%-14s","L2rho","L2u");
	if (DIM >= 2)
		index += sprintf(header_spec+index,"%-14s","L2v");
	if (DIM >= 3)
		index += sprintf(header_spec+index,"%-14s","L2w");
	sprintf(header_spec+index,"%-14s%-14s","L2p","L2s");

	return header_spec;
}

static const char* compute_header_spec_euler_entropy ( )
{
	static char header_spec[STRLEN_MAX];
	sprintf(header_spec,"%-14s","L2s");
	return header_spec;
}
