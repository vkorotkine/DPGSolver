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

#include "macros.h"
#include "definitions_bc.h"

#include "def_templates_face_solver_opg.h"
#include "def_templates_face_solver.h"
#include "def_templates_volume_solver.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_compute_face_rlhs.h"
#include "def_templates_penalty_opg.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for the inverse mass matrix of the input face.
 *  \return See brief. */
static const struct const_Matrix_T* constructor_inverse_mass_face_T
	(const struct Solver_Face_T*const s_face, ///< Standard.
	 const struct const_Matrix_T*const mass   ///< Face mass matrix. Input if available, otherwise pass `NULL`.
	);

/// \brief Set the function pointers for \ref OPG_Solver_Face_T::constructor_rlhs_penalty.
static void set_function_pointers_penalty_T
	(struct OPG_Solver_Face_T*const opg_s_face ///< Standard.
	 );

// Interface functions ********************************************************************************************** //

void constructor_derived_OPG_Solver_Face_T (struct Face* face_ptr, const struct Simulation* sim)
{
	UNUSED(sim);
	struct Solver_Face_T* s_face         = (struct Solver_Face_T*) face_ptr;
	struct OPG_Solver_Face_T* opg_s_face = (struct OPG_Solver_Face_T*) face_ptr;

	opg_s_face->m_inv = constructor_inverse_mass_face_T(s_face,NULL); // destructed

	set_function_pointers_penalty_T(opg_s_face);
}

void destructor_derived_OPG_Solver_Face_T (struct Face* face_ptr)
{
	struct OPG_Solver_Face_T* opg_s_face = (struct OPG_Solver_Face_T*) face_ptr;

	destructor_const_Matrix_T(opg_s_face->m_inv);
}

const struct Operator* get_operator__cv0_vt_fc_T
	(const int side_index, const struct OPG_Solver_Face_T*const opg_s_face)
{
	const struct Face*const face             = (struct Face*) opg_s_face;
	const struct Solver_Face_T*const s_face  = (struct Solver_Face_T*) opg_s_face;
	const struct Volume*const vol            = face->neigh_info[side_index].volume;
	const struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) vol;
	const struct OPG_Solver_Element*const e  = (struct OPG_Solver_Element*) vol->element;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = face->neigh_info[side_index].ind_href;
	const int p_v = s_vol->p_ref,
	          p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	return get_Multiarray_Operator(e->cv0_vt_fc[curved],(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
}

struct Multiarray_Operator get_operator__cv1_vt_fc_T
	(const int side_index, const struct OPG_Solver_Face_T*const opg_s_face)
{
	const struct Face*const face             = (struct Face*) opg_s_face;
	const struct Solver_Face_T*const s_face  = (struct Solver_Face_T*) opg_s_face;
	const struct Volume*const vol            = face->neigh_info[side_index].volume;
	const struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) vol;
	const struct OPG_Solver_Element*const e  = (struct OPG_Solver_Element*) vol->element;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = face->neigh_info[side_index].ind_href;
	const int p_v = s_vol->p_ref,
	          p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	return set_MO_from_MO(e->cv1_vt_fc[curved],1,(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Constructor for the mass matrix of the input face.
 *  \return See brief. */
static const struct const_Matrix_T* constructor_mass_face_T
	(const struct Solver_Face_T*const s_face ///< Standard.
	);

/// \brief Set the function pointers for \ref OPG_Solver_Face_T::constructor_rlhs_penalty for boundary faces.
static void set_function_pointers_penalty_boundary_T
	(struct OPG_Solver_Face_T*const opg_s_face ///< Standard.
	 );

static const struct const_Matrix_T* constructor_inverse_mass_face_T
	(const struct Solver_Face_T*const s_face, const struct const_Matrix_T*const mass)
{
	const struct const_Matrix_T* m_inv = NULL;
	if (mass) {
		m_inv = constructor_inverse_const_Matrix_T(mass); // returned
	} else {
		const struct const_Matrix_T*const mass = constructor_mass_face_T(s_face); // destructed
		m_inv = constructor_inverse_const_Matrix_T(mass); // returned
		destructor_const_Matrix_T(mass);
	}
	return m_inv;
}

static void set_function_pointers_penalty_T (struct OPG_Solver_Face_T*const opg_s_face)
{
	const struct Face*const face = (struct Face*) opg_s_face;
	if (!face->boundary) {
		opg_s_face->constructor_rlhs_penalty[0] = constructor_rlhs_f_test_penalty_unsupported_T;
		opg_s_face->constructor_rlhs_penalty[1] = constructor_rlhs_f_test_penalty_unsupported_T;
	} else {
		set_function_pointers_penalty_boundary_T(opg_s_face);
	}
}

// Level 1 ********************************************************************************************************** //

/** \brief Set the function pointers for \ref OPG_Solver_Face_T::constructor_rlhs_penalty for boundary faces when
 * solving the linear advection equation. */
static void set_function_pointers_penalty_boundary_advection_T
	(struct OPG_Solver_Face_T*const opg_s_face ///< Standard.
	 );

static const struct const_Matrix_T* constructor_mass_face_T (const struct Solver_Face_T*const s_face)
{
	const struct Operator*const cv0_ff_fc  = get_operator__cv0_ff_fc_T(s_face);
	const struct const_Vector_R*const w_fc = get_operator__w_fc__s_e_T(s_face);
	const struct const_Vector_T jac_det_fc = interpret_const_Multiarray_as_Vector_T(s_face->jacobian_det_fc);

	const struct const_Vector_T*const wJ_fc = constructor_dot_mult_const_Vector_T_RT(1.0,w_fc,&jac_det_fc,1); // dest.

	const struct const_Matrix_R*const m_l = cv0_ff_fc->op_std;
	const struct const_Matrix_T*const m_r = constructor_mm_diag_const_Matrix_R_T(1.0,m_l,wJ_fc,'L',false); // destructed
	destructor_const_Vector_T(wJ_fc);

	const struct const_Matrix_T*const mass = constructor_mm_RT_const_Matrix_T('T','N',1.0,m_l,m_r,'R'); // returned
	destructor_const_Matrix_T(m_r);

	return mass;
}

static void set_function_pointers_penalty_boundary_T (struct OPG_Solver_Face_T*const opg_s_face)
{
	switch (get_set_pde_index(NULL)) {
	case PDE_ADVECTION:
		set_function_pointers_penalty_boundary_advection_T(opg_s_face);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",get_set_pde_index(NULL));
		break;
	}
}

// Level 2 ********************************************************************************************************** //

static void set_function_pointers_penalty_boundary_advection_T (struct OPG_Solver_Face_T*const opg_s_face)
{
	const struct Face* face = (struct Face*) opg_s_face;
	const int bc = face->bc % BC_STEP_SC;

	opg_s_face->constructor_rlhs_penalty[0] = constructor_rhs_f_test_penalty_do_nothing_T;
	switch (bc) {
	case BC_UPWIND:      case BC_UPWIND_ALT1: case BC_UPWIND_ALT2:
	case BC_UPWIND_ALT3: case BC_UPWIND_ALT4: case BC_UPWIND_ALT5:
		opg_s_face->constructor_rlhs_penalty[1] = constructor_rlhs_f_test_penalty_advection_upwind_T;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",face->bc);
		break;
	}
}

#include "undef_templates_face_solver_opg.h"
#include "undef_templates_face_solver.h"
#include "undef_templates_volume_solver.h"

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_compute_face_rlhs.h"
#include "undef_templates_penalty_opg.h"
