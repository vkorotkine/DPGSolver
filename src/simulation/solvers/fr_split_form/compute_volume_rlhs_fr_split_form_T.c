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
#include "definitions_intrusive.h"


#include "def_templates_compute_volume_rlhs_fr_split_form.h"

#include "def_templates_volume_solver.h"
#include "def_templates_volume_solver_fr_split_form.h"

#include "def_templates_compute_volume_rlhs.h"
#include "def_templates_flux.h"
#include "def_templates_test_case.h"
#include "def_templates_operators.h"

#include "def_templates_matrix.h"
#include "element_operators.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"
// Static function declarations ************************************************************************************* //

/// \brief Container for solver-related parameters.
struct S_Params_T {
	struct S_Params_Volume_Structor_T spvs; ///< \ref S_Params_Volume_Structor_T.

	compute_rlhs_v_fptr_T compute_rlhs; ///< Pointer to the appropriate function.
};

/** \brief Set the parameters of \ref S_Params_T.
 *  \return A statically allocated \ref S_Params_T container. */
static struct S_Params_T set_s_params_T
	(const struct Simulation* sim ///< \ref Simulation.
	);
static void compute_rhs_v_frsf_like_T
	(const struct Flux_Ref_T*const flux_r, struct Solver_Volume_T*const s_vol,
	 struct Solver_Storage_Implicit*const ssi);

/* static const struct const_Matrix_d* constructor_mass_frsf_T (const struct Solver_Volume_T* s_vol); */
// Interface functions ********************************************************************************************** //

void compute_volume_rlhs_fr_split_form_T
	(const struct Simulation* sim, struct Solver_Storage_Implicit* ssi, struct Intrusive_List* volumes)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER_FRSF);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_FRSF);

	struct S_Params_T s_params = set_s_params_T(sim);
	struct Flux_Input_T* flux_i = constructor_Flux_Input_T(sim); // destructed

	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) curr;

		struct Flux_Ref_T* flux_r = constructor_Flux_Ref_vol_T(&s_params.spvs,flux_i,s_vol); // destructed

		// Compute the rhs (and optionally the lhs) terms.
		s_params.compute_rlhs(flux_r,s_vol,ssi);
		destructor_Flux_Ref_T(flux_r);
	}
	destructor_Flux_Input_T(flux_i);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct S_Params_T set_s_params_T (const struct Simulation* sim)
{
	struct S_Params_T s_params;

	set_S_Params_Volume_Structor_T(&s_params.spvs,sim);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	switch (test_case->solver_method_curr) {
	case 'e':
		s_params.compute_rlhs = compute_rhs_v_frsf_like_T;
		break;
	case 'i':
	default:
		EXIT_ERROR("Unsupported: %c (type_rc: %d)\n",test_case->solver_method_curr,TYPE_RC);
		break;
	}

	return s_params;
}
static void compute_rhs_v_frsf_like_T
	(const struct Flux_Ref_T*const flux_r, struct Solver_Volume_T*const s_vol,
	 struct Solver_Storage_Implicit*const ssi)
{
	UNUSED(flux_r);
	UNUSED(ssi);

	struct Multiarray_d* sol_coef = s_vol->sol_coef;

	ptrdiff_t size = compute_size(sol_coef->order,sol_coef->extents);
	struct Vector_d* u  = constructor_empty_Vector_d(size);
	struct Vector_d* u2 = constructor_empty_Vector_d(size);
	for (int i = 0; i < size; ++i) {
		u->data[i] = sol_coef->data[i];
		u2->data[i] = 0.5*sol_coef->data[i]*sol_coef->data[i];
	}
	/* print_Vector_d(u); */
	struct Volume* vol = (struct Volume*) s_vol;
	struct FRSF_Solver_Element* e = (struct FRSF_Solver_Element*) vol->element;
	const struct const_Matrix_d* D = e->D_frsf;
	/* print_const_Matrix_d(e->D_frsf); */
	const struct const_Vector_d* Du2 = constructor_mv_const_Vector_d('N',-2.0/3.0,D,(struct const_Vector_d*)u2);
	const struct const_Vector_d* Du  = constructor_mv_const_Vector_d('N',1.0,D,(struct const_Vector_d*)u);
	const struct const_Vector_d* uDu = constructor_dot_mult_const_Vector_d(-1.0/3.0,(struct const_Vector_d*)u,Du,1);

	/* print_Vector_d(u); */
	/* print_Vector_d(u2); */
	/* print_const_Vector_d(Du2); */
	/* print_const_Vector_d(Du); */
	/* print_const_Vector_d(uDu); */

	const struct const_Vector_d* sum = constructor_sum_Vectors_const_Vector_d(1.0,Du2,1.0,uDu);
	/* print_const_Vector_d(sum); */

	for (int i = 0; i < sum->ext_0; ++i)
		s_vol->rhs->data[i] = sum->data[i];
	/* print_Multiarray_d(s_vol->rhs); */
	struct FRSF_Solver_Volume_T* frsf_s_vol = (struct FRSF_Solver_Volume_T*) s_vol;
	/* frsf_s_vol->m     = constructor_mass_frsf_T(s_vol); */
	/* printf("mass elem %lf\n",frsf_s_vol->m->data[4]); */
	/* print_Multiarray_T(s_vol->rhs); */
	/* print_const_Matrix_T(frsf_s_vol->m);//confirmed mass matrix perf */
	for (int i = 0; i < sum->ext_0; ++i){
		/* double k=(sum->ext_0)*i; */
		/* int j; */
		/* j = (int)(k); */
		int j = 9*i;
		s_vol->rhs->data[i] = frsf_s_vol->m->data[j]*s_vol->rhs->data[i];//multiply rhs by M so scale M_inv later

		/* printf("minv val %lf\n",frsf_s_vol->m_inv->data[j]); */
	}
	/* print_Multiarray_d(s_vol->rhs); */
	/* EXIT; */
}
/* static const struct const_Matrix_d* constructor_mass_frsf_T (const struct Solver_Volume_T* s_vol) */
/* { */
/* 	const struct Operator*const cv0_vs_vc  = get_operator__cv0_vs_vc(s_vol); */
/* 	const struct const_Vector_d*const w_vc = get_operator__w_vc__s_e(s_vol); */

/* 	const struct const_Matrix_d*const m_l = cv0_vs_vc->op_std; */
/* 	const struct const_Matrix_d*const m_r = constructor_mm_diag_const_Matrix_d(1.0,m_l,w_vc,'L',false); // destructed */

/* 	const struct const_Matrix_d*const mass = constructor_mm_const_Matrix_d('T','N',1.0,m_l,m_r,'R'); // returned */
/* 	destructor_const_Matrix_d(m_r); */

/* 	return mass; */
/* } */
#include "undef_templates_compute_volume_rlhs_fr_split_form.h"

#include "undef_templates_volume_solver.h"
#include "undef_templates_volume_solver_fr_split_form.h"

#include "undef_templates_compute_volume_rlhs.h"
#include "undef_templates_flux.h"
#include "undef_templates_test_case.h"
#include "undef_templates_operators.h"
