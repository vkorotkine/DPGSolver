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
#include <string.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_test_case.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"
#include "def_templates_volume_solver.h"
#include "def_templates_volume_solver_dg.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/// Container holding flags for which members of \ref DG_Solver_Volume_T are needed.
struct Needed_Members {
	bool sol_coef_p, ///< Flag for \ref DG_Solver_Volume_T::sol_coef_p.
	     m_inv,      ///< Flag for \ref DG_Solver_Volume_T::m_inv.
	     m;          ///< Flag for \ref DG_Solver_Volume_T::m.
};

/** \brief Return a statically allocated \ref Needed_Members container with values set.
 *  \return See brief. */
static struct Needed_Members set_needed_members
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for the mass matrix of the input volume.
 *  \return See brief. */
static const struct const_Matrix_R* constructor_mass
	(const struct Solver_Volume_T* s_vol ///< \ref Solver_Volume_T.
	);

/** \brief Constructor for the inverse mass matrix of the input volume.
 *  \return See brief. */
static const struct const_Matrix_R* constructor_inverse_mass
	(const struct DG_Solver_Volume_T* dg_s_vol ///< \ref DG_Solver_Volume_T.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_DG_Solver_Volume_T (struct Volume* volume_ptr, const struct Simulation* sim)
{
	struct Needed_Members needed_members = set_needed_members(sim);

	struct Solver_Volume_T* s_vol       = (struct Solver_Volume_T*) volume_ptr;
	struct DG_Solver_Volume_T* dg_s_vol = (struct DG_Solver_Volume_T*) volume_ptr;

	const int order = s_vol->sol_coef->order;
	ptrdiff_t* extents = s_vol->sol_coef->extents;

	dg_s_vol->rhs        = constructor_empty_Multiarray_T('C',order,extents); // destructed
	dg_s_vol->sol_coef_p =
		( needed_members.sol_coef_p ? constructor_empty_Multiarray_T('C',order,extents) : NULL ); // destructed

	dg_s_vol->m     = ( needed_members.m     ? constructor_mass(s_vol) : NULL );            // destructed
	dg_s_vol->m_inv = ( needed_members.m_inv ? constructor_inverse_mass(dg_s_vol) : NULL ); // destructed

	const struct Test_Case_T*const test_case = (struct Test_Case_T*) sim->test_case_rc->tc;
	if (test_case->has_2nd_order) {
		const int order = s_vol->grad_coef->order;
		ptrdiff_t* extents = s_vol->grad_coef->extents;
		dg_s_vol->grad_coef_v = constructor_zero_Multiarray_T('C',order,extents); // destructed

		if (test_case->solver_method_curr == 'i') {
			for (int i = 0; i < DIM; ++i)
				dg_s_vol->d_g_coef_v__d_s_coef[i] = constructor_empty_const_Matrix_R('R',0,0); // destructed
		} else {
			assert(test_case->solver_method_curr == 'e');
		}
	}
}

void destructor_derived_DG_Solver_Volume_T (struct Volume* volume_ptr)
{
	struct DG_Solver_Volume_T* dg_s_vol = (struct DG_Solver_Volume_T*) volume_ptr;

	destructor_Multiarray_T(dg_s_vol->rhs);
	destructor_conditional_Multiarray_T(dg_s_vol->sol_coef_p);
	destructor_conditional_const_Matrix_R(dg_s_vol->m_inv);
	destructor_conditional_const_Matrix_R(dg_s_vol->m);
	destructor_conditional_Multiarray_T(dg_s_vol->grad_coef_v);
	for (int i = 0; i < DIM; ++i)
		destructor_conditional_const_Matrix_R(dg_s_vol->d_g_coef_v__d_s_coef[i]);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Needed_Members set_needed_members (const struct Simulation* sim)
{
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	struct Needed_Members needed_members =
		{ .sol_coef_p = false,
		  .m          = false,
		  .m_inv      = false, };

	switch (test_case->solver_proc) {
	case SOLVER_E: // fallthrough
	case SOLVER_EI:
		if (!sim->collocated)
			needed_members.m_inv = true;
		switch (test_case->solver_type_e) {
		case SOLVER_E_SSP_RK_33: // fallthrough
		case SOLVER_E_LS_RK_54:
			needed_members.sol_coef_p = true;
			break;
		case SOLVER_E_EULER:
			break; // Do nothing
		default:
			EXIT_ERROR("Unsupported: %d\n",test_case->solver_type_e);
			break;
		}
		break;
	case SOLVER_I:
		// Do nothing
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->solver_proc);
		break;
	}

	if (test_case->has_2nd_order) {
		switch (sim->method) {
		case METHOD_DG:
			if (!sim->collocated)
				needed_members.m_inv = true;
			break;
		case METHOD_DPG:
			break; // Do nothing.
		default:
			EXIT_ERROR("Unsupported: %d",sim->method);
			break;
		}
	}

	if (test_case->lhs_terms == LHS_CFL_RAMPING)
		needed_members.m = true;

	return needed_members;
}

static const struct const_Matrix_R* constructor_mass (const struct Solver_Volume_T* s_vol)
{
	struct Volume* vol               = (struct Volume*) s_vol;
	const struct Solver_Element* s_e = (struct Solver_Element*) vol->element;

	const int p      = s_vol->p_ref,
	          curved = vol->curved;
	const struct Operator* cv0_vs_vc = get_Multiarray_Operator(s_e->cv0_vs_vc[curved],(ptrdiff_t[]){0,0,p,p});
	const struct const_Vector_R* w_vc = get_operator__w_vc__s_e_T(s_vol);

	const struct const_Vector_R jacobian_det_vc = interpret_const_Multiarray_as_Vector_R(s_vol->jacobian_det_vc);
	const struct const_Vector_R* wJ_vc = constructor_dot_mult_const_Vector_R(1.0,w_vc,&jacobian_det_vc,1); // destructed

	const struct const_Matrix_R* m_l = cv0_vs_vc->op_std;
	const struct const_Matrix_R* m_r = constructor_mm_diag_const_Matrix_R(1.0,m_l,wJ_vc,'L',false); // destructed
	destructor_const_Vector_R(wJ_vc);

	const struct const_Matrix_R* mass = constructor_mm_const_Matrix_R('T','N',1.0,m_l,m_r,'R'); // returned
	destructor_const_Matrix_R(m_r);

	return mass;
}

static const struct const_Matrix_R* constructor_inverse_mass (const struct DG_Solver_Volume_T* dg_s_vol)
{
	const struct const_Matrix_R* m_inv = NULL;
	if (dg_s_vol->m) {
		m_inv = constructor_inverse_const_Matrix_R(dg_s_vol->m); // returned
	} else {
		const struct const_Matrix_R* m = constructor_mass((struct Solver_Volume_T*)dg_s_vol); // destructed
		m_inv = constructor_inverse_const_Matrix_R(m); // returned
		destructor_const_Matrix_R(m);
	}
	return m_inv;
}
