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
#include "definitions_numerical_flux.h"


#include "def_templates_compute_face_rlhs_opg.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_face_solver.h"
#include "def_templates_face_solver_opg.h"
#include "def_templates_volume_solver.h"

#include "def_templates_boundary.h"
#include "def_templates_compute_face_rlhs.h"
#include "def_templates_flux.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_operators.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for solver-related parameters.
struct S_Params_T {
	compute_rlhs_opg_f_fptr_T compute_rlhs; ///< Pointer to the appropriate function.
};

/// \brief Container for numerical flux related parameters.
struct Num_Flux_T {
	const struct const_Multiarray_T* n_dot_nf; ///< Unit normal dotted with the numerical flux.
};

/** \brief Set the parameters of \ref S_Params_T.
 *  \return A statically allocated \ref S_Params_T container. */
static struct S_Params_T set_s_params_T
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Constructor for a \ref Numerical_Flux_T container holding the correct values for the OPG method.
 *  \return See brief.
 *
 *  For internal faces, the normal numerical flux values are simply interpolated from the appropriate coefficients of
 *  \ref Solver_Face_T. For boundary faces, the normal numerical flux values are computed as for the DG scheme.
 */
static struct Numerical_Flux_T* constructor_Numerical_Flux_OPG_T
	(struct Numerical_Flux_Input_T*const num_flux_i, ///< Standard.
	 const struct Solver_Face_T*const s_face,        ///< Standard.
	 const struct Simulation*const sim               ///< Standard.
	);

/** \brief Constructor for a \ref Flux_T container holding the correct values for the OPG method.
 *  \return On boundary faces, returns the constructed \ref Flux_T; otherwise returns NULL. */
static struct Flux_T* constructor_Flux_OPG_T
	(struct Flux_Input_T*const flux_i,       ///< Standard.
	 const struct Solver_Face_T*const s_face ///< Standard.
	);

/// \brief Scale \ref Numerical_Flux_T::nnf by the face Jacobian (i.e. only the explicit term).
static void scale_by_Jacobian_e_T
	(struct Numerical_Flux_T*const num_flux, ///< See brief.
	 const struct Solver_Face_T*const s_face ///< See brief.
	);

/** \brief Return the jump of the solution test functions, [[v]], at the face cubature nodes.
 *  \return See brief.
 *
 *  The jump is defined by:
 *  - internal faces: [[v]] := v_l - v_r;
 *  - boundary faces: [[v]] := v_l;
 */
static const struct const_Multiarray_T* constructor_jump_test_s_fc_T
	(const struct OPG_Solver_Face_T*const opg_s_face ///< Standard.
		);

// Interface functions ********************************************************************************************** //

void compute_face_rlhs_opg_T
	(const struct Simulation*const sim, struct Solver_Storage_Implicit*const ssi, struct Intrusive_List*const faces)
{
	if (get_set_has_1st_2nd_order(NULL)[1])
		EXIT_ERROR("Add support/Ensure that all is working as expected.\n");

	assert(sim->elements->name == IL_ELEMENT_SOLVER_OPG);
	assert(sim->faces->name    == IL_FACE_SOLVER_OPG);
	assert(sim->volumes->name  == IL_VOLUME_SOLVER_OPG);

	struct S_Params_T s_params = set_s_params_T(sim);

	struct Test_Case_T*const test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	assert(test_case->solver_method_curr == 'i');

	struct Numerical_Flux_Input_T*const num_flux_i = constructor_Numerical_Flux_Input_T(sim); // destructed
	struct Flux_Input_T*const flux_i               = constructor_Flux_Input_T(sim); // destructed

	for (struct Intrusive_Link* curr = faces->first; curr; curr = curr->next) {
		struct Solver_Face_T*const s_face = (struct Solver_Face_T*) curr;

		struct Numerical_Flux_T*const num_flux = constructor_Numerical_Flux_OPG_T(num_flux_i,s_face,sim); // dest.
		struct Flux_T*const flux               = constructor_Flux_OPG_T(flux_i,s_face); // destructed

		scale_by_Jacobian_e_T(num_flux,s_face);
		s_params.compute_rlhs(flux,num_flux,s_face,ssi);
		destructor_Numerical_Flux_T(num_flux);
		destructor_conditional_Flux_T(flux);
	}
	destructor_Numerical_Flux_Input_T(num_flux_i);
	destructor_Flux_Input_T(flux_i);
}

void update_coef_nf_f_opg_T (const struct Simulation*const sim)
{
	/** The L2 projection is current being used to project the test to the normal flux coefficients. For face
	 *  collocated schemes when the normal flux polynomial and test function polynomial degrees are equal, the
	 *  projection operator reduces to identity. */
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Face*const face                  = (struct Face*) curr;
		if (!UPDATE_NF_BOUNDARY && face->boundary)
			continue;

		const struct Solver_Face_T*const s_face         = (struct Solver_Face_T*) curr;
		const struct OPG_Solver_Face_T*const opg_s_face = (struct OPG_Solver_Face_T*) curr;

		const struct Operator*const cv0_ff_fc = get_operator__cv0_ff_fc(0,s_face);
		const struct const_Matrix_T*const op_proj_L2 =
			constructor_mm_TR_const_Matrix_T('N','T',1.0,opg_s_face->m_inv,cv0_ff_fc->op_std,'R'); // dest.

		const struct const_Vector_R*const w_fc  = get_operator__w_fc__s_e(s_face);
		const struct const_Vector_T j_det_fc    = interpret_const_Multiarray_as_Vector_T(s_face->jacobian_det_fc);
		const struct const_Vector_T*const wJ_fc = constructor_dot_mult_const_Vector_T_RT(1.0,w_fc,&j_det_fc,1); // d.

		const struct const_Multiarray_T*const jump_test_s_fc = constructor_jump_test_s_fc_T(opg_s_face); // dest.
		scale_Multiarray_by_Vector_T('L',1.0,(struct Multiarray_T*)jump_test_s_fc,wJ_fc,false);
		destructor_const_Vector_T(wJ_fc);

		mm_NNC_Multiarray_T(1.0,0.0,op_proj_L2,jump_test_s_fc,s_face->nf_coef);
		destructor_const_Matrix_T(op_proj_L2);
		destructor_const_Multiarray_T(jump_test_s_fc);
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Version of \ref compute_rlhs_opg_f_fptr_T used to call \ref compute_rhs_f_dg_like_T.
 *  \return See brief. */
static void compute_rhs_f_opg_dg_like_T
	(const struct Flux_T*const flux,               ///< Standard.
	 const struct Numerical_Flux_T*const num_flux, ///< Standard.
	 struct Solver_Face_T*const s_face,            ///< Standard.
	 struct Solver_Storage_Implicit*const ssi      ///< Standard.
		);

/** \brief Construct the data members of the \ref Numerical_Flux_Input_T container which are specific to the face under
 *         consideration for the opg scheme. */
static void constructor_Numerical_Flux_Input_data_opg_T
	(struct Numerical_Flux_Input_T*const num_flux_i, ///< Standard.
	 const struct Solver_Face_T*const s_face,        ///< Standard.
	 const struct Simulation*const sim               ///< Standard.
	);

static struct S_Params_T set_s_params_T (const struct Simulation*const sim)
{
	struct S_Params_T s_params;

	struct Test_Case_T*const test_case = (struct Test_Case_T*)sim->test_case_rc->tc;

	switch (test_case->solver_method_curr) {
/// \todo Uncomment the preprocessor directives once the complex version of the function is added.
//#if TYPE_RC == TYPE_COMPLEX
	case 'e':
		s_params.compute_rlhs = compute_rhs_f_opg_dg_like_T;
		break;
//#elif TYPE_RC == TYPE_REAL
#if TYPE_RC == TYPE_REAL
	case 'i':
		if (test_case->has_1st_order && !test_case->has_2nd_order) {
			s_params.compute_rlhs = compute_rlhs_1;
		} else if (!test_case->has_1st_order && test_case->has_2nd_order) {
			EXIT_ADD_SUPPORT;
		} else if (test_case->has_1st_order && test_case->has_2nd_order) {
			EXIT_ADD_SUPPORT;
		} else {
			EXIT_ERROR("Unsupported: %d %d\n",test_case->has_1st_order,test_case->has_2nd_order);
		}
		break;
#endif
	default:
		EXIT_ERROR("Unsupported: %c (type_rc: %d)\n",test_case->solver_method_curr,TYPE_RC);
		break;
	}

	return s_params;
}

static struct Numerical_Flux_T* constructor_Numerical_Flux_OPG_T
	(struct Numerical_Flux_Input_T*const num_flux_i, const struct Solver_Face_T*const s_face,
	 const struct Simulation*const sim)
{
	struct Numerical_Flux* num_flux = NULL;

	const struct Face*const face = (struct Face*) s_face;
	if (!face->boundary) {
		const struct Operator*const cv0_ff_fc = get_operator__cv0_ff_fc_T(0,s_face);
		const struct const_Multiarray_T*const nf_coef = (struct const_Multiarray_T*) s_face->nf_coef;

		num_flux = calloc(1,sizeof *num_flux); // returned
		num_flux->nnf =
			constructor_mm_NN1_Operator_const_Multiarray_T(cv0_ff_fc,nf_coef,'C','d',nf_coef->order,NULL); // d.
	} else {
		constructor_Numerical_Flux_Input_data_opg_T(num_flux_i,s_face,sim); // destructed
		num_flux = constructor_Numerical_Flux_T(num_flux_i); // returned
		destructor_Numerical_Flux_Input_data_T(num_flux_i);
	}
	return num_flux;
}

static struct Flux_T* constructor_Flux_OPG_T (struct Flux_Input_T*const flux_i, const struct Solver_Face_T*const s_face)
{
	const struct Face*const face = (struct Face*) s_face;
	if (!face->boundary)
		return NULL;

	constructor_Flux_Input_data_f_T(0,flux_i,s_face);
	struct Flux_T*const flux = constructor_Flux_T(flux_i); // returned
	destructor_Flux_Input_data_f_T(flux_i);

	return flux;
}

static void scale_by_Jacobian_e_T
	(struct Numerical_Flux_T*const num_flux, const struct Solver_Face_T*const s_face)
{
	const struct const_Vector_R jacobian_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	scale_Multiarray_T_by_Vector_R('L',1.0,(struct Multiarray_T*)num_flux->nnf,&jacobian_det_fc,false);
}

static const struct const_Multiarray_T* constructor_jump_test_s_fc_T (const struct OPG_Solver_Face_T*const opg_s_face)
{
	const struct Face*const face            = (struct Face*) opg_s_face;
	const struct Solver_Face_T*const s_face = (struct Solver_Face_T*) opg_s_face;
	const struct const_Matrix_R* cv0_vt_fc[2] = { get_operator__cv0_vt_fc_T(0,opg_s_face)->op_std, NULL, };
	const struct Solver_Volume_T* s_vol[2]    = { (struct Solver_Volume_T*) face->neigh_info[0].volume, NULL, };

	const struct const_Multiarray_T*const test_s_coef = (struct const_Multiarray_T*) s_vol[0]->test_s_coef;
	struct Multiarray_T*const jump_test_s = constructor_mm_NN1C_Multiarray_T(cv0_vt_fc[0],test_s_coef); // returned

	if (!face->boundary) {
		cv0_vt_fc[1] = constructor_copy_const_Matrix_R(get_operator__cv0_vt_fc_T(1,opg_s_face)->op_std); // dest.
		permute_Matrix_R_fc((struct Matrix_R*)cv0_vt_fc[1],'R',0,s_face);
		s_vol[1] = (struct Solver_Volume_T*) face->neigh_info[1].volume;

		const struct const_Multiarray_T*const test_s_coef = (struct const_Multiarray_T*) s_vol[1]->test_s_coef;
		const struct const_Multiarray_T*const jump_test_s_R =
			constructor_mm_NN1C_const_Multiarray_T(cv0_vt_fc[1],test_s_coef); // destructed
		destructor_const_Matrix_R(cv0_vt_fc[1]);

		add_in_place_Multiarray_T(-1.0,jump_test_s,jump_test_s_R);
		destructor_const_Multiarray_T(jump_test_s_R);
	}

	return (struct const_Multiarray_T*) jump_test_s;
}

// Level 1 ********************************************************************************************************** //

static void constructor_Numerical_Flux_Input_data_opg_T
	(struct Numerical_Flux_Input_T*const num_flux_i, const struct Solver_Face_T*const s_face,
	 const struct Simulation*const sim)
{
	if (get_set_has_1st_2nd_order(NULL)[1])
		EXIT_ADD_SUPPORT;
	constructor_Numerical_Flux_Input_data_T(num_flux_i,s_face,sim); // destructed
}

static void compute_rhs_f_opg_dg_like_T
	(const struct Flux_T*const flux, const struct Numerical_Flux_T*const num_flux, struct Solver_Face_T*const s_face,
	 struct Solver_Storage_Implicit*const ssi)
{
	UNUSED(flux);
	compute_rhs_f_dg_like_T(num_flux,s_face,ssi);
}

#include "undef_templates_compute_face_rlhs_opg.h"

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_face_solver.h"
#include "undef_templates_face_solver_opg.h"
#include "undef_templates_volume_solver.h"

#include "undef_templates_boundary.h"
#include "undef_templates_compute_face_rlhs.h"
#include "undef_templates_flux.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_operators.h"
#include "undef_templates_test_case.h"
