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

#include "compute_face_rlhs_opg.h"

#include "definitions_test_case.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "face_solver_opg.h"
#include "element_solver_opg.h"
#include "volume.h"
#include "volume_solver_opg.h"

#include "computational_elements.h"
#include "compute_rlhs.h"
#include "compute_face_rlhs.h"
#include "multiarray_operator.h"
#include "numerical_flux.h"
#include "operator.h"
#include "simulation.h"
#include "solve.h"
#include "solve_opg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Version of \ref compute_rlhs_f_fptr_T computing the rhs and lhs terms for 1st order equations only.
static void compute_rlhs_1
	(const struct Numerical_Flux*const num_flux, ///< See brief.
	 struct Solver_Face*const s_face,            ///< See brief.
	 struct Solver_Storage_Implicit*const ssi    ///< See brief.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "compute_face_rlhs_opg_T.c"

void update_coef_nf_f_opg (const struct Simulation*const sim)
{
	/** The L2 projection is current being used to project the test to the normal flux coefficients. For face
	 *  collocated schemes when the normal flux polynomial and test function polynomial degrees are equal, the
	 *  projection operator reduces to identity. */
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		/* 1. Compute test jump on the face (special case for boundary faces).
		 * 2. Compute the L2 projection operator.
		 * 3. Update nf_coef.
		 */
UNUSED(curr); EXIT_ADD_SUPPORT;
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref compute_rlhs_f_fptr_T computing the lhs terms for 1st order equations only.
static void compute_lhs_1
	(struct Solver_Face*const s_face,         ///< See brief.
	 struct Solver_Storage_Implicit*const ssi ///< See brief.
	);

static void compute_rlhs_1
	(const struct Numerical_Flux*const num_flux, struct Solver_Face*const s_face,
	 struct Solver_Storage_Implicit*const ssi)
{
	compute_rhs_f_dg_like(num_flux,s_face,ssi);
	compute_lhs_1(s_face,ssi);
}

// Level 1 ********************************************************************************************************** //

/// \brief Container for operators needed for the assembly of LHS terms for the OPG scheme.
struct Lhs_Operators_OPG {
	const struct const_Vector_d* wJ_fc; ///< Face cubature weights "dot-multiplied" by the Jacobian determinant.

	/** 'c'oefficient to 'v'alue operators from the 'v'olume 't'est basis to 'f'ace 'c'ubature nodes.
	 *
	 *  The indices are used to denote the following operators:
	 *  - [0]: ll operator ('l'eft  basis -> 'l'eft  nodes);
	 *  - [1]: rl operator ('r'ight basis -> 'l'eft  nodes; permutation of node ordering);
	 */
	const struct const_Matrix_d* cv_vt_fc[2];
};

/** \brief Constructor for a \ref Lhs_Operators_OPG container.
 *  \return See brief. */
static const struct Lhs_Operators_OPG* constructor_Lhs_Operators_OPG
	(const struct Solver_Face*const s_face ///< Standard.
	);

/// \brief Destructor for a \ref Lhs_Operators_OPG container.
static void destructor_Lhs_Operators_OPG
	(const struct Solver_Face*const s_face, ///< Standard.
	 const struct Lhs_Operators_OPG* ops    ///< Standard.
	);

/// \brief Finalize the 1st order lhs term contribution from the \ref Face for the opg scheme.
static void finalize_lhs_1_f_opg
	(const int side_index[2],                  /**< The indices of the affectee, affector, respectively. See the
	                                            *   comments in \ref compute_face_rlhs_dg.h for the convention. */
	 const struct Lhs_Operators_OPG*const ops, ///< Standard.
	 const struct Solver_Face*const s_face,    ///< Standard.
	 struct Solver_Storage_Implicit*const ssi  ///< Standard.
	);

static void compute_lhs_1 (struct Solver_Face*const s_face, struct Solver_Storage_Implicit*const ssi)
{
	/// See \ref compute_face_rlhs_dg.h for the `lhs_**` notation.
	const struct Face*const face = (struct Face*) s_face;

	const struct Lhs_Operators_OPG*const ops = constructor_Lhs_Operators_OPG(s_face); // destructed

	finalize_lhs_1_f_opg((int[]){0,0},ops,s_face,ssi); // lhs_ll
	if (!face->boundary) {
		finalize_lhs_1_f_opg((int[]){0,1},ops,s_face,ssi); // lhs_lr
		finalize_lhs_1_f_opg((int[]){1,0},ops,s_face,ssi); // lhs_rl
		finalize_lhs_1_f_opg((int[]){1,1},ops,s_face,ssi); // lhs_rr
	}
	destructor_Lhs_Operators_OPG(s_face,ops);
}

// Level 2 ********************************************************************************************************** //

/** \brief Get the pointer to the appropriate \ref OPG_Solver_Element::cv0_vt_fc operator.
 *  \return See brief. */
static const struct Operator* get_operator__cv0_vt_fc
	(const int side_index,                 ///< The index of the side of the face under consideration.
	 const struct Solver_Face*const s_face ///< The current \ref Face.
	);

static const struct Lhs_Operators_OPG* constructor_Lhs_Operators_OPG (const struct Solver_Face*const s_face)
{
	struct Lhs_Operators_OPG*const ops = calloc(1,sizeof *ops); // free

	const struct const_Vector_d*const w_fc  = get_operator__w_fc__s_e(s_face);
	const struct const_Vector_d j_det_fc    = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	ops->wJ_fc = constructor_dot_mult_const_Vector_d(1.0,w_fc,&j_det_fc,1); // destructed

	const struct Operator*const cv0_vt_fc = get_operator__cv0_vt_fc(0,s_face);
	ops->cv_vt_fc[0] = cv0_vt_fc->op_std;

	const struct Face*const face = (struct Face*) s_face;
	if (!face->boundary) {
		ops->cv_vt_fc[1] = constructor_copy_const_Matrix_d(get_operator__cv0_vt_fc(1,s_face)->op_std); // dest.
		permute_Matrix_d_fc((struct Matrix_d*)ops->cv_vt_fc[1],'R',0,s_face);
	}
	return ops;
}

static void destructor_Lhs_Operators_OPG (const struct Solver_Face*const s_face, const struct Lhs_Operators_OPG* ops)
{
	destructor_const_Vector_d(ops->wJ_fc);

	const struct Face*const face = (struct Face*) s_face;
	if (!face->boundary)
		destructor_const_Matrix_d(ops->cv_vt_fc[1]);
	free((void*)ops);
}

static void finalize_lhs_1_f_opg
	(const int side_index[2], const struct Lhs_Operators_OPG*const ops, const struct Solver_Face*const s_face,
	 struct Solver_Storage_Implicit*const ssi)
{
	const struct Face*const face = (struct Face*) s_face;
	const struct OPG_Solver_Volume*const opg_s_vol[2] = { (struct OPG_Solver_Volume*) face->neigh_info[0].volume,
	                                                      (struct OPG_Solver_Volume*) face->neigh_info[1].volume, };

	const struct const_Matrix_d*const lhs_r =
		constructor_mm_diag_const_Matrix_d_d(1.0,ops->cv_vt_fc[side_index[1]],ops->wJ_fc,'L',false); // destructed.

	const double scale = ( side_index[0] == side_index[1] ? 1.0 : -1.0 );
	const struct const_Matrix_d*const lhs =
		constructor_mm_const_Matrix_d('T','N',scale,ops->cv_vt_fc[side_index[0]],lhs_r,'R'); // destructed
	destructor_const_Matrix_d(lhs_r);

	const int*const n_vr_eq = get_set_n_var_eq(NULL);
	const int n_vr    = n_vr_eq[0],
	          n_eq    = n_vr_eq[1];
	/** \warning It is possible that a change may be required when systems of equations are used. Currently, there is
	 *           a "default coupling" between the face terms between each equations and variables. */
	 assert(n_vr == 1 && n_eq == 1); // Ensure that all is working properly when removed.

	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		set_petsc_Mat_row_col_opg(ssi,opg_s_vol[side_index[0]],eq,opg_s_vol[side_index[1]],vr);
		add_to_petsc_Mat(ssi,lhs);
	}}
}

// Level 3 ********************************************************************************************************** //

static const struct Operator* get_operator__cv0_vt_fc (const int side_index, const struct Solver_Face*const s_face)
{
	const struct Face*const face            = (struct Face*) s_face;
	const struct Volume*const vol           = face->neigh_info[side_index].volume;
	const struct Solver_Volume*const s_vol  = (struct Solver_Volume*) vol;
	const struct OPG_Solver_Element*const e = (struct OPG_Solver_Element*) vol->element;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = face->neigh_info[side_index].ind_href;
	const int p_v = s_vol->p_ref,
	          p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	return get_Multiarray_Operator(e->cv0_vt_fc[curved],(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
}
