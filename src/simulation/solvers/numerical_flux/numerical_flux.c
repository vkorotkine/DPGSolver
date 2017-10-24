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

#include "numerical_flux.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

#include "multiarray.h"
#include "vector.h"

#include "element_solver_dg.h"
#include "face.h"
#include "face_solver.h"
#include "volume.h"
#include "volume_solver.h"

#include "const_cast.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve_dg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

struct Numerical_Flux_Input* constructor_Numerical_Flux_Input (const struct Simulation* sim)
{
	struct Numerical_Flux_Input* num_flux_i = calloc(1,sizeof *num_flux_i); // returned

	struct Test_Case* test_case = sim->test_case;
	const_cast_i(&num_flux_i->d,sim->d);
	const_cast_i(&num_flux_i->n_eq,test_case->n_eq);
	const_cast_i(&num_flux_i->n_var,test_case->n_var);

	const_cast_b(&num_flux_i->has_1st_order,test_case->has_1st_order);
	const_cast_b(&num_flux_i->has_2nd_order,test_case->has_2nd_order);

	num_flux_i->compute_Numerical_Flux = test_case->compute_Numerical_Flux;
	switch (test_case->solver_method_curr) {
	case 'e':
		num_flux_i->compute_member = test_case->flux_comp_mem_e;
		num_flux_i->compute_Numerical_Flux_1st = test_case->compute_Numerical_Flux_e[0];
		num_flux_i->compute_Numerical_Flux_2nd = test_case->compute_Numerical_Flux_e[1];
		break;
	case 'i':
		num_flux_i->compute_member = test_case->flux_comp_mem_i;
		num_flux_i->compute_Numerical_Flux_1st = test_case->compute_Numerical_Flux_i[0];
		num_flux_i->compute_Numerical_Flux_2nd = test_case->compute_Numerical_Flux_i[1];
		break;
	default:
		EXIT_ERROR("Unsupported: %c.\n",test_case->solver_method_curr);
		break;
	}

	return num_flux_i;
}

void destructor_Numerical_Flux_Input (struct Numerical_Flux_Input* num_flux_i)
{
	free(num_flux_i);
}

struct Numerical_Flux* constructor_Numerical_Flux (const struct Numerical_Flux_Input* num_flux_i)
{
	assert(((num_flux_i->neigh_info[0].s != NULL && num_flux_i->neigh_info[0].s->layout == 'C') &&
	        (num_flux_i->neigh_info[1].s != NULL && num_flux_i->neigh_info[1].s->layout == 'C')) ||
	       ((num_flux_i->neigh_info[0].g != NULL && num_flux_i->neigh_info[0].g->layout == 'C') &&
	        (num_flux_i->neigh_info[1].g != NULL && num_flux_i->neigh_info[1].g->layout == 'C')));

	const bool* c_m = num_flux_i->compute_member;
	assert(c_m[0] || c_m[1] || c_m[2]);

	const int d     = num_flux_i->d,
	          n_eq  = num_flux_i->n_eq,
		    n_var = num_flux_i->n_var;
	const ptrdiff_t n_n = ( num_flux_i->neigh_info[0].s != NULL ? num_flux_i->neigh_info[0].s->extents[0] :
	                                                              num_flux_i->neigh_info[0].g->extents[0] );

	struct mutable_Numerical_Flux* num_flux = calloc(1,sizeof *num_flux); // returned

	num_flux->nnf = (c_m[0] ? constructor_zero_Multiarray_d('C',2,(ptrdiff_t[]){n_n,n_eq}) : NULL);
	for (int i = 0; i < 2; ++i) {
		struct m_Neigh_Info_NF n_i = num_flux->neigh_info[i];
		n_i.dnnf_ds = (c_m[1] ? constructor_zero_Multiarray_d('C',3,(ptrdiff_t[]){n_n,n_eq,n_var})   : NULL);
		n_i.dnnf_dg = (c_m[2] ? constructor_zero_Multiarray_d('C',4,(ptrdiff_t[]){n_n,n_eq,n_var,d}) : NULL);
	}

	num_flux_i->compute_Numerical_Flux(num_flux_i,num_flux);

	return (struct Numerical_Flux*) num_flux;
}

void destructor_Numerical_Flux (struct Numerical_Flux* num_flux)
{
	if (num_flux->nnf != NULL)
		destructor_const_Multiarray_d(num_flux->nnf);

	for (int i = 0; i < 2; ++i) {
		struct Neigh_Info_NF n_i = num_flux->neigh_info[i];
		if (n_i.dnnf_ds != NULL)
			destructor_const_Multiarray_d(n_i.dnnf_ds);
		if (n_i.dnnf_dg != NULL)
			destructor_const_Multiarray_d(n_i.dnnf_dg);
	}
}

void compute_Numerical_Flux_1 (const struct Numerical_Flux_Input* num_flux_i, struct mutable_Numerical_Flux* num_flux)
{
	num_flux_i->compute_Numerical_Flux_1st(num_flux_i,num_flux);
}

void compute_Numerical_Flux_12 (const struct Numerical_Flux_Input* num_flux_i, struct mutable_Numerical_Flux* num_flux)
{
	num_flux_i->compute_Numerical_Flux_1st(num_flux_i,num_flux);
	num_flux_i->compute_Numerical_Flux_2nd(num_flux_i,num_flux);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Version of \ref constructor_sg_fc_fptr interpolating the solution from the neighbouring volume to the face.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_s_fc_interp
	(const struct Face* face,      ///< Defined for \ref constructor_sg_fc_fptr.
	 const struct Simulation* sim, ///< Defined for \ref constructor_sg_fc_fptr.
	 const int side_index          ///< The index of the side of the face under consideration.
	);

const struct const_Multiarray_d* constructor_s_l_fcl_interp (const struct Face* face, const struct Simulation* sim)
{
	return constructor_s_fc_interp(face,sim,0);
}

const struct const_Multiarray_d* constructor_s_r_fcl_interp (const struct Face* face, const struct Simulation* sim)
{
	const int side_index = 1;

	struct Multiarray_d* sol_r_fcr = (struct Multiarray_d*) constructor_s_fc_interp(face,sim,side_index);
	permute_Multiarray_d_fc(sol_r_fcr,'R',side_index,(struct Solver_Face*)face);

	return (const struct const_Multiarray_d*) sol_r_fcr;
}

const struct const_Multiarray_d* constructor_sg_fc_null (const struct Face* face, const struct Simulation* sim)
{
	UNUSED(face);
	UNUSED(sim);
	return NULL;
}

// Level 1 ********************************************************************************************************** //

static const struct const_Multiarray_d* constructor_s_fc_interp
	(const struct Face* face, const struct Simulation* sim, const int side_index)
{
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	struct Solver_Face* s_face     = (struct Solver_Face*) face;
	struct Volume* volume          = face->neigh_info[side_index].volume;
	struct Solver_Volume* s_volume = (struct Solver_Volume*) volume;

	const struct DG_Solver_Element* e = (const struct DG_Solver_Element*) volume->element;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = face->neigh_info[side_index].ind_href;
	const int p_v = s_volume->p_ref,
	          p_f = s_face->p_ref;

	const struct Operator* cv0_vs_fc = ( (s_face->cub_type == 's')
		? get_Multiarray_Operator(e->cv0_vs_fcs,(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v})
		: get_Multiarray_Operator(e->cv0_vs_fcc,(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v}) );

	const struct const_Multiarray_d* s_coef = (const struct const_Multiarray_d*) s_volume->sol_coef;

	return constructor_mm_NN1_Operator_const_Multiarray_d(cv0_vs_fc,s_coef,'C',op_format,s_coef->order,NULL);
}
