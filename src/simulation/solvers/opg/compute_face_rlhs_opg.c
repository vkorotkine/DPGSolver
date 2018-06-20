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
#include "volume_solver_opg.h"

#include "computational_elements.h"
#include "compute_rlhs.h"
#include "compute_face_rlhs.h"
#include "compute_volume_rlhs.h"
#include "compute_volume_rlhs_opg.h"
#include "multiarray_operator.h"
#include "numerical_flux.h"
#include "operator.h"
#include "simulation.h"
#include "solve.h"
#include "solve_opg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** Flag for whether the normal flux on the boundary faces should be updated. At the time of this implementation,
 *  \ref Solver_Face_T::nf_coef was not stored on the boundary faces as boundary normal fluxes were computed using the
 *  numerical flux; it may consequently be required to initialize nf_coef on boundary faces if enabling this option. */
#define UPDATE_NF_BOUNDARY false

/// \brief Version of \ref compute_rlhs_opg_f_fptr_T computing the rhs and lhs terms for 1st order equations only.
static void compute_rlhs_1
	(const struct Flux_Ref*const flux_r,         ///< See brief.
	 const struct Numerical_Flux*const num_flux, ///< See brief.
	 struct Solver_Face*const s_face,            ///< See brief.
	 struct Solver_Storage_Implicit*const ssi    ///< See brief.
		);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "compute_face_rlhs_opg_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "compute_face_rlhs_opg_T.c"
#include "undef_templates_type.h"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Version of \ref compute_rlhs_f_fptr_T computing the lhs terms for 1st order equations only for internal
 *         faces. */
static void compute_lhs_1_i
	(struct OPG_Solver_Face*const opg_s_face, ///< See brief.
	 struct Solver_Storage_Implicit*const ssi ///< See brief.
	);

/** \brief Version of \ref compute_rlhs_f_fptr_T computing the lhs terms for 1st order equations only for boundary
 *         faces. */
static void compute_lhs_1_b
	(const struct Flux_Ref*const flux_r,         ///< See brief.
	 const struct Numerical_Flux*const num_flux, ///< See brief.
	 struct Solver_Face*const s_face,            ///< See brief.
	 struct Solver_Storage_Implicit*const ssi    ///< See brief.
		);

static void compute_rlhs_1
	(const struct Flux_Ref*const flux_r, const struct Numerical_Flux*const num_flux, struct Solver_Face*const s_face,
	 struct Solver_Storage_Implicit*const ssi)
{
	compute_rhs_f_dg_like(num_flux,s_face,ssi);

	struct OPG_Solver_Face*const opg_s_face = (struct OPG_Solver_Face*) s_face;
	const struct Face*const face = (struct Face*) s_face;
	if (!face->boundary)
		compute_lhs_1_i(opg_s_face,ssi);
	else
		compute_lhs_1_b(flux_r,num_flux,s_face,ssi);
	opg_s_face->constructor_rlhs_penalty(flux_r,num_flux,s_face,ssi);
}

// Level 1 ********************************************************************************************************** //

/// \brief Finalize the 1st order lhs term contribution from the \ref Face for the opg scheme.
static void finalize_lhs_1_f_opg
	(const int side_index[2],                  /**< The indices of the affectee, affector, respectively. See the
	                                            *   comments in \ref compute_face_rlhs_dg.h for the convention. */
	 const struct Lhs_Operators_OPG*const ops, ///< Standard.
	 const struct Solver_Face*const s_face,    ///< Standard.
	 struct Solver_Storage_Implicit*const ssi  ///< Standard.
	);

/** \brief Constructor for the lhs term arising from the OPG scheme for boundary faces.
 *  \return See brief.
 *
 *  This term is nearly identical to that for the DG scheme with the modification being that the linearization is with
 *  respect to the test function coefficients and not the solution coefficients.
 */
static const struct const_Matrix_d* constructor_lhs_f_1_b
	(const struct Flux_Ref*const flux_ref,         ///< Standard.
	 const struct Numerical_Flux*const num_flux, ///< Standard.
	 const struct Solver_Face*const s_face       ///< Standard.
	);

static void compute_lhs_1_i (struct OPG_Solver_Face*const opg_s_face, struct Solver_Storage_Implicit*const ssi)
{
	/// See \ref compute_face_rlhs_dg.h for the `lhs_**` notation.
	const struct Face*const face = (struct Face*) opg_s_face;
	assert(!face->boundary);

	const struct Solver_Face*const s_face = (struct Solver_Face*) opg_s_face;

	const struct Lhs_Operators_OPG*const ops = constructor_Lhs_Operators_OPG(opg_s_face); // destructed

	finalize_lhs_1_f_opg((int[]){0,0},ops,s_face,ssi); // lhs_ll
	finalize_lhs_1_f_opg((int[]){0,1},ops,s_face,ssi); // lhs_lr
	finalize_lhs_1_f_opg((int[]){1,0},ops,s_face,ssi); // lhs_rl
	finalize_lhs_1_f_opg((int[]){1,1},ops,s_face,ssi); // lhs_rr

	destructor_Lhs_Operators_OPG(s_face,ops);
}

static void compute_lhs_1_b
	(const struct Flux_Ref*const flux_r, const struct Numerical_Flux*const num_flux, struct Solver_Face*const s_face,
	 struct Solver_Storage_Implicit*const ssi)
{
	const struct Face*const face = (struct Face*) s_face;
	assert(face->boundary);

	const struct const_Matrix_d*const lhs = constructor_lhs_f_1_b(flux_r,num_flux,s_face); // destructed

	const struct OPG_Solver_Volume*const opg_s_vol = (struct OPG_Solver_Volume*) face->neigh_info[0].volume;
	set_petsc_Mat_row_col_opg(ssi,opg_s_vol,0,opg_s_vol,0);
	add_to_petsc_Mat(ssi,lhs);

	destructor_const_Matrix_d(lhs);
}

// Level 2 ********************************************************************************************************** //

/** \brief Constructor for the 'l'eft term of the lhs contribution from the boundary faces.
 *  \return The operator corresponding to test_s' frac{dnnf_ds}. */
static const struct const_Matrix_d* constructor_lhs_f_1_b_l
	(const struct Numerical_Flux*const num_flux, ///< Standard.
	 const struct Solver_Face*const s_face       ///< Standard.
	);

/** \brief Constructor for the 'r'ight term of the lhs contribution from the boundary faces.
 *  \return The operator corresponding to
 *          frac{d s}{d test_s_coef} = cv0_vs_fc M^{-1} cv0_vt_vc' (-df_ds' (dot) cv1_vt_vc). */
static const struct const_Matrix_d* constructor_lhs_f_1_b_r
	(const struct Flux_Ref*const flux_r,   ///< Standard.
	 const struct Solver_Face*const s_face ///< Standard.
	);

static void finalize_lhs_1_f_opg
	(const int side_index[2], const struct Lhs_Operators_OPG*const ops, const struct Solver_Face*const s_face,
	 struct Solver_Storage_Implicit*const ssi)
{
	const struct Face*const face = (struct Face*) s_face;
	const struct OPG_Solver_Volume*const opg_s_vol[2] = { (struct OPG_Solver_Volume*) face->neigh_info[0].volume,
	                                                      (struct OPG_Solver_Volume*) face->neigh_info[1].volume, };

	/** The linearization is determined according to the relation between \ref Solver_Face_T::nf_coef and
	 *  \ref Solver_Volume_T::test_s_coef. */
	const struct const_Matrix_d*const lhs_l_p1 =
		constructor_mm_diag_const_Matrix_d_d(1.0,ops->cv0_vt_fc[side_index[0]],ops->wJ_fc,'L',false); // dest.
	const struct const_Matrix_d*const cv0_ff_fc = get_operator__cv0_ff_fc(s_face)->op_std;
	const struct const_Matrix_d*const lhs_l_p2 = constructor_mm_const_Matrix_d('T','N',1.0,lhs_l_p1,cv0_ff_fc,'R'); // d.
	destructor_const_Matrix_d(lhs_l_p1);

	const struct const_Matrix_d*const lhs_l = constructor_mm_const_Matrix_d('N','N',1.0,lhs_l_p2,ops->proj_L2_l,'R'); // d.
	destructor_const_Matrix_d(lhs_l_p2);

#if 0 // Only enable for face collocated parameters (p_t_p = 0).
	const struct const_Multiarray_d*const normals_fc = s_face->normals_fc;
	const ptrdiff_t n_fc = normals_fc->extents[0];
	const double b_adv[] = { 0.5, 0.5, 0.0, };
	struct Vector_d* b_dot_n2_V = constructor_zero_Vector_d(n_fc); // destructed.
	double* b_dot_n2 = b_dot_n2_V->data;
	for (int i = 0; i < n_fc; ++i) {
		for (int j = 0; j < DIM; ++j)
			b_dot_n2[i] += b_adv[j]*get_row_const_Multiarray_d(i,normals_fc)[j];
//		b_dot_n2[i] = sqrt(b_dot_n2[i]*b_dot_n2[i]);
	}
	scale_Matrix_by_Vector_d('R',1.0,(struct Matrix_d*)lhs_l,(struct const_Vector_d*)b_dot_n2_V,false);
	destructor_Vector_d(b_dot_n2_V);
#endif

	const double scale = ( side_index[0] == side_index[1] ? -1.0 : 1.0 );
	const struct const_Matrix_d*const lhs =
		constructor_mm_const_Matrix_d('N','N',scale,lhs_l,ops->cv0_vt_fc[side_index[1]],'R'); // destructed
	destructor_const_Matrix_d(lhs_l);

	const int*const n_vr_eq = get_set_n_var_eq(NULL);
	const int n_vr = n_vr_eq[0];
	const int n_eq = n_vr_eq[1];
	/** \warning It is possible that a change may be required when systems of equations are used. Currently, there is
	 *           a "default coupling" between the face terms between each equations and variables.
	 *
	 *  From a few of the DPG papers, it seems that an additional sign(n (dot) df/du) scaling may need to be added.
	 *  While this would introduce the coupling between the equations, it may also destroy the symmetry for non
	 *  scalar PDEs. Entropy variables to fix this or simply don't assume symmetric? THINK.
	 */
	 assert(n_vr == 1 && n_eq == 1); // Ensure that all is working properly when removed.

	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		if (eq != vr)
			continue;
		set_petsc_Mat_row_col_opg(ssi,opg_s_vol[side_index[0]],eq,opg_s_vol[side_index[1]],vr);
		add_to_petsc_Mat(ssi,lhs);
	}}
}

static const struct const_Matrix_d* constructor_lhs_f_1_b
	(const struct Flux_Ref*const flux_r, const struct Numerical_Flux*const num_flux,
	 const struct Solver_Face*const s_face)
{
	const struct const_Matrix_d*const lhs_l  = constructor_lhs_f_1_b_l(num_flux,s_face); // destructed
	const struct const_Matrix_d*const lhs_r  = constructor_lhs_f_1_b_r(flux_r,s_face);   // destructed

	const struct const_Matrix_d*const lhs = constructor_mm_const_Matrix_d('N','N',1.0,lhs_l,lhs_r,'R'); // returned
	destructor_const_Matrix_d(lhs_l);
	destructor_const_Matrix_d(lhs_r);

	return lhs;
}

// Level 3 ********************************************************************************************************** //

static const struct const_Matrix_d* constructor_lhs_f_1_b_l
	(const struct Numerical_Flux*const num_flux, const struct Solver_Face*const s_face)
{
	const struct Face*const face = (struct Face*) s_face;
	assert(face->boundary);

	const int side_index[2] = { 0, 0, };
	struct Matrix_d*const lhs_l = constructor_lhs_f_1(side_index,num_flux,s_face); // returned

	return (struct const_Matrix_d*) lhs_l;
}

static const struct const_Matrix_d* constructor_lhs_f_1_b_r
	(const struct Flux_Ref*const flux_r, const struct Solver_Face*const s_face)
{
	const struct Face*const face = (struct Face*) s_face;
	assert(face->boundary);

	const struct OPG_Solver_Volume*const opg_s_vol = (struct OPG_Solver_Volume*) face->neigh_info[0].volume;

	return constructor_operator__test_s_coef_to_sol_coef_d(flux_r,opg_s_vol);
}
