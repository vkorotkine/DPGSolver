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

#include "def_templates_face_solver_opg.h"
#include "def_templates_face_solver.h"
#include "def_templates_volume_solver.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_compute_face_rlhs.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for the inverse mass matrix of the input face.
 *  \return See brief. */
static const struct const_Matrix_R* constructor_inverse_mass_face_T
	(const struct Solver_Face_T*const s_face, ///< Standard.
	 const struct const_Matrix_R*const mass   ///< Face mass matrix. Input if available, otherwise pass `NULL`.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_OPG_Solver_Face_T (struct Face* face_ptr, const struct Simulation* sim)
{
	UNUSED(sim);
	struct Solver_Face_T* s_face         = (struct Solver_Face_T*) face_ptr;
	struct OPG_Solver_Face_T* opg_s_face = (struct OPG_Solver_Face_T*) face_ptr;

	opg_s_face->m_inv = constructor_inverse_mass_face_T(s_face,NULL); // destructed
}

void destructor_derived_OPG_Solver_Face_T (struct Face* face_ptr)
{
	struct OPG_Solver_Face_T* opg_s_face = (struct OPG_Solver_Face_T*) face_ptr;

	destructor_const_Matrix_R(opg_s_face->m_inv);
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
static const struct const_Matrix_R* constructor_mass_face_T
	(const struct Solver_Face_T*const s_face ///< Standard.
	);

static const struct const_Matrix_R* constructor_inverse_mass_face_T
	(const struct Solver_Face_T*const s_face, const struct const_Matrix_R*const mass)
{
	const struct const_Matrix_R* m_inv = NULL;
	if (mass) {
		m_inv = constructor_inverse_const_Matrix_R(mass); // returned
	} else {
		const struct const_Matrix_R*const mass = constructor_mass_face_T(s_face); // destructed
		m_inv = constructor_inverse_const_Matrix_R(mass); // returned
		destructor_const_Matrix_R(mass);
	}
	return m_inv;
}

// Level 1 ********************************************************************************************************** //

static const struct const_Matrix_R* constructor_mass_face_T (const struct Solver_Face_T*const s_face)
{
	const struct Operator*const cv0_ff_fc  = get_operator__cv0_ff_fc_T(0,s_face);
	const struct const_Vector_R*const w_fc = get_operator__w_fc__s_e_T(s_face);
	const struct const_Vector_R jac_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);

	const struct const_Vector_R*const wJ_fc = constructor_dot_mult_const_Vector_R(1.0,w_fc,&jac_det_fc,1); // dest.

	const struct const_Matrix_R*const m_l = cv0_ff_fc->op_std;
	const struct const_Matrix_R*const m_r = constructor_mm_diag_const_Matrix_R(1.0,m_l,wJ_fc,'L',false); // destructed
	destructor_const_Vector_R(wJ_fc);

	const struct const_Matrix_R*const mass = constructor_mm_const_Matrix_R('T','N',1.0,m_l,m_r,'R'); // returned
	destructor_const_Matrix_R(m_r);

	return mass;
}

#include "undef_templates_face_solver_opg.h"
#include "undef_templates_face_solver.h"
#include "undef_templates_volume_solver.h"

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_compute_face_rlhs.h"
