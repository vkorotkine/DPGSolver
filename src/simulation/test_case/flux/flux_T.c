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

#include "macros.h"
#include "definitions_core.h"


#include "def_templates_flux.h"

#include "def_templates_multiarray.h"

#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

struct Flux_Input_T* constructor_Flux_Input_T (const struct Simulation* sim)
{
	struct Flux_Input_T* flux_i = calloc(1,sizeof *flux_i); // returned

	const_cast_c1(&flux_i->input_path,sim->input_path);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const_cast_i(&flux_i->n_eq,test_case->n_eq);
	const_cast_i(&flux_i->n_var,test_case->n_var);

	const_cast_b(&flux_i->has_1st_order,test_case->has_1st_order);
	const_cast_b(&flux_i->has_2nd_order,test_case->has_2nd_order);

	flux_i->compute_Flux = test_case->compute_Flux;
	flux_i->compute_Flux_1st = test_case->compute_Flux_iv[0];
	flux_i->compute_Flux_2nd = test_case->compute_Flux_iv[1];
	switch (test_case->solver_method_curr) {
	case 'e':
		flux_i->compute_member = test_case->flux_comp_mem_e;
		break;
	case 'i':
		flux_i->compute_member = test_case->flux_comp_mem_i;
		break;
	default:
		EXIT_ERROR("Unsupported: %c.\n",test_case->solver_method_curr);
		break;
	}

	return flux_i;
}

void destructor_Flux_Input_T (struct Flux_Input_T* flux_i)
{
	free(flux_i);
}

struct Flux_T* constructor_Flux_T (const struct Flux_Input_T* flux_i)
{
	assert((flux_i->s != NULL && flux_i->s->layout == 'C') ||
	       (flux_i->g != NULL && flux_i->g->layout == 'C'));

	const bool* compute_member = flux_i->compute_member;
	assert(compute_member[0]);

	const int n_eq = flux_i->n_eq,
	          n_vr = flux_i->n_var;
	const ptrdiff_t n_n = ( flux_i->s != NULL ? flux_i->s->extents[0] : flux_i->g->extents[0] );

	struct mutable_Flux_T* flux = calloc(1,sizeof *flux); // returned

	flux->f     = (compute_member[0] ?
		constructor_zero_Multiarray_T('C',3,(ptrdiff_t[]){n_n,DIM,n_eq})           : NULL); // destructed
	flux->df_ds = (compute_member[1] ?
		constructor_zero_Multiarray_T('C',4,(ptrdiff_t[]){n_n,DIM,n_eq,n_vr})      : NULL); // destructed
	flux->df_dg = (compute_member[2] ?
		constructor_zero_Multiarray_T('C',5,(ptrdiff_t[]){n_n,DIM,n_eq,n_vr,DIM})  : NULL); // destructed
	flux->d2f_ds2 = (compute_member[3] ?
		constructor_zero_Multiarray_T('C',5,(ptrdiff_t[]){n_n,DIM,n_eq,n_vr,n_vr}) : NULL); // destructed

	flux_i->compute_Flux(flux_i,flux);

	return (struct Flux_T*) flux;
}

void destructor_Flux_T (struct Flux_T* flux)
{
	if (flux->f != NULL)
		destructor_const_Multiarray_T(flux->f);
	if (flux->df_ds != NULL)
		destructor_const_Multiarray_T(flux->df_ds);
	if (flux->df_dg != NULL)
		destructor_const_Multiarray_T(flux->df_dg);
	if (flux->d2f_ds2 != NULL)
		destructor_const_Multiarray_T(flux->d2f_ds2);

	free(flux);
}

void compute_Flux_1_T (const struct Flux_Input_T* flux_i, struct mutable_Flux_T* flux)
{
	flux_i->compute_Flux_1st(flux_i,flux);
}

void compute_Flux_12_T (const struct Flux_Input_T* flux_i, struct mutable_Flux_T* flux)
{
	flux_i->compute_Flux_1st(flux_i,flux);
	flux_i->compute_Flux_2nd(flux_i,flux);
}

void increment_pointers_T (const int n_ptr, const Type**const ptrs)
{
	for (int i = 0; i < n_ptr; ++i)
		++ptrs[i];
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
