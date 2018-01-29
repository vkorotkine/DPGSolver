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
#include "definitions_test_case.h"

#include "def_templates_compute_grad_coef_dg.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_face_solver.h"
#include "def_templates_face_solver_dg.h"
#include "def_templates_volume_solver.h"
#include "def_templates_volume_solver_dg.h"

#include "def_templates_boundary.h"
#include "def_templates_compute_face_rlhs.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_solve.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Compute \ref DG_Solver_Volume_T::grad_coef_v for the list of volumes.
static void compute_grad_coef_volumes
	(const struct Simulation*const sim,  ///< \ref Simulation.
	 struct Intrusive_List*const volumes ///< The list of volumes.
	);

/// \brief Compute \ref DG_Solver_Face_T::grad_coef_f for the list of faces.
static void compute_grad_coef_faces
	(const struct Simulation*const sim, ///< \ref Simulation.
	 struct Intrusive_List*const faces  ///< The list of faces.
	);

// Interface functions ********************************************************************************************** //

void compute_grad_coef_dg_T
	(const struct Simulation*const sim, struct Intrusive_List*const volumes, struct Intrusive_List*const faces)
{
	const struct Test_Case_T*const test_case = (struct Test_Case_T*) sim->test_case_rc->tc;
	if (!test_case->has_2nd_order)
		return;

	compute_grad_coef_volumes(sim,volumes);
	compute_grad_coef_faces(sim,faces);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Get the appropriate sub-range of the \ref DG_Solver_Element::cv1_vs_vc operators.
 *  \return See brief. */
static struct Multiarray_Operator get_operator__cv1_vs_vc
	(const struct DG_Solver_Volume_T*const dg_s_vol ///< \ref DG_Solver_Volume_T.
	);

/** \brief Constructor for the partial physical gradient operator from the reference operator and metric terms.
 *  \return See brief.
 *
 *  A partial operator is returned in the sense that the inverse Jacobian determinant contribution is omitted.
 */
static const struct const_Matrix_R* constructor_grad_xyz_p
	(const int dir,                                   ///< The direction. Options: 0 (x), 1 (y), 2 (z).
	 const struct Multiarray_Operator*const grad_rst, ///< The reference gradient operators.
	 const struct const_Multiarray_R*const metrics    ///< The metric terms.
	);

/// \brief Exit with an error if the numerical solution is not central for the current \ref Test_Cast_T::ind_num_flux.
static void assert_numerical_solution_is_central
	(const int ind_num_flux_2nd ///< The 2nd \ref Test_Cast_T::ind_num_flux.
	);

/** \brief Constructor for the difference of numerical solution and interpolated solution from the left volume at the
 *         face cubature nodes.
 *  \return See brief. */
static const struct const_Matrix_T* constructor_diff_s_num_s
	(const int ind_num_flux_2nd,                    ///< Defined for \ref compute_d_g_coef_f__d_s_coef_i.
	 const struct DG_Solver_Face_T*const dg_s_face, ///< \ref DG_Solver_Face_T.
	 const struct Simulation*const sim              ///< \ref Simulation.
	);

/** \brief Constructor for \ref const_Matrix_T\*s holding the Jacobian determinant, normal and solution difference
 *         terms.
 *  \return See brief. */
static const struct const_Matrix_T*const* constructor_jdet_n_diff_fc
	(const struct const_Matrix_R*const jdet_n_fc,   ///< Jacobian det. dot multiplied with normal vector terms.
	 const struct const_Matrix_T*const diff_s_num_s ///< Return from \ref constructor_diff_s_num_s.
	);

/// \brief Destructor for \ref constructor_jdet_n_diff_fc.
static void destructor_jdet_n_diff_fc
	(const struct const_Matrix_T*const* jdet_n_diff_fc ///< Standard.
	);

/// \brief Computes the \ref DG_Solver_Face_T::Neigh_Info_DG::g_coef_f term associated with the current side.
static void compute_g_coef_f_i
	(const int side_index,                                   ///< Defined for \ref get_sol_scale.
	 const struct const_Matrix_T*const*const jdet_n_diff_fc, /**< Jacobian det., normals and solution difference at
	                                                          *   face cubature nodes of the current side. */
	 struct DG_Solver_Face_T*const dg_s_face,                ///< \ref DG_Solver_Face_T.
	 const bool collocated                                   ///< \ref Simulation::collocated.
	);

/** \brief Computes the \ref DG_Solver_Face_T::Neigh_Info_DG::d_g_coef_f__d_s_coef term associated with the current side
 *         for 'i'nternal faces.
 *
 *  The `side_index` input determines the \ref DG_Solver_Face::neigh_info.
 */
static void compute_d_g_coef_f__d_s_coef_i
	(const int side_index,                        ///< Defined for \ref get_sol_scale.
	 const int ind_num_flux_2nd,                  ///< Defined for \ref get_sol_scale.
	 const struct const_Matrix_R*const jdet_n_fc, ///< Jacobian determinant dotted with normals at fc nodes.
	 struct DG_Solver_Face_T*const dg_s_face,     ///< \ref DG_Solver_Face_T.
	 const bool collocated                        ///< \ref Simulation::collocated.
	);

/** \brief Computes the \ref DG_Solver_Face_T::Neigh_Info_DG::g_coef_f term associated with the current side using the
 *         associated linearization terms.
 *
 *  \warning It is assumed that the gradients are linear wrt to the solution in this function.
 */
static void compute_g_coef_f_i_using_lin
	(const int side_index,                   ///< The index of the side under consideration.
	 struct DG_Solver_Face_T*const dg_s_face ///< \ref DG_Solver_Face_T.
	);

/** \brief Computes gradient-related members from \ref DG_Solver_Face_T for boundary faces.
 *
 *  The following members are computed:
 *  - DG_Solver_Face_T::Neigh_Info_DG::grad_coef_f;
 *  - DG_Solver_Face_T::Neigh_Info_DG::d_g_coef_f__d_s_coef (side_index = 0, sol_index = 0 only).
 */
static void compute_g_coef_related_boundary
	(const int ind_num_flux_2nd,                  ///< Defined for \ref compute_d_g_coef_f__d_s_coef_i.
	 const struct const_Matrix_R*const jdet_n_fc, ///< Defined for \ref compute_d_g_coef_f__d_s_coef_i.
	 struct DG_Solver_Face_T*const dg_s_face,     ///< Defined for \ref compute_d_g_coef_f__d_s_coef_i.
	 const struct Simulation*const sim            ///< \ref Simulation.
	);

/// \brief Add the \ref DG_Solver_Face_T::Neigh_Info_DG::g_coef_f contributions to \ref Solver_Volume_T::grad_coef.
static void add_face_grad_coef_f_to_volumes
	(const struct DG_Solver_Face_T*const dg_s_face ///< \ref DG_Solver_Face_T.
	);

static void compute_grad_coef_volumes (const struct Simulation*const sim, struct Intrusive_List*const volumes)
{
	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
		const struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) curr;
		struct DG_Solver_Volume_T*const dg_s_vol = (struct DG_Solver_Volume_T*) curr;

		for (int d = 0; d < DIM; ++d) {
			const struct Multiarray_Operator cv1_vs_vc = get_operator__cv1_vs_vc(dg_s_vol);
			const struct const_Matrix_R*const grad_xyz =
				constructor_grad_xyz_p(d,&cv1_vs_vc,s_vol->metrics_vc); // destructed/keep

			const struct const_Matrix_R* d_g_coef_v__d_s_coef = NULL;
			if (!sim->collocated) {
				const struct Operator*const tw0_vt_vc = get_operator__tw0_vt_vc_T(s_vol);
				const struct const_Matrix_R*const ibp2_v =
					constructor_mm_const_Matrix_R('N','N',1.0,tw0_vt_vc->op_std,grad_xyz,'R'); // destructed
				d_g_coef_v__d_s_coef =
					constructor_mm_const_Matrix_R('N','N',1.0,dg_s_vol->m_inv,ibp2_v,'R'); // keep
				destructor_const_Matrix_R(ibp2_v);
			} else {
				// Multiplication by the cubature weights is **not** performed here as multiplication by the
				// inverse cubature weights is **not** performed when multiplying by the inverse mass matrix.
				const struct const_Vector_R jacobian_det_vc =
					interpret_const_Multiarray_as_Vector_R(s_vol->jacobian_det_vc);
				scale_Matrix_R_by_Vector_R('L',1.0,(struct Matrix_R*)grad_xyz,&jacobian_det_vc,true);
				d_g_coef_v__d_s_coef = grad_xyz;
			}

			struct Multiarray_T grad_coef_v =
				interpret_Multiarray_as_slice_T(dg_s_vol->grad_coef_v,2,(ptrdiff_t[]){d});
			mm_NNC_Multiarray_T(1.0,0.0,d_g_coef_v__d_s_coef,
			                    (struct const_Multiarray_T*)s_vol->sol_coef,&grad_coef_v);

			const struct Test_Case_T*const test_case = (struct Test_Case_T*) sim->test_case_rc->tc;
			if (test_case->solver_method_curr == 'i') {
				destructor_const_Matrix_R(dg_s_vol->d_g_coef_v__d_s_coef[d]);
				dg_s_vol->d_g_coef_v__d_s_coef[d] = d_g_coef_v__d_s_coef;
				if (!sim->collocated)
					destructor_const_Matrix_R(grad_xyz);
			} else {
				assert(test_case->solver_method_curr == 'e');
				destructor_const_Matrix_R(grad_xyz);
			}
		}
		copy_into_Multiarray_T(s_vol->grad_coef,(struct const_Multiarray_T*)dg_s_vol->grad_coef_v);
	}
}

static void compute_grad_coef_faces (const struct Simulation*const sim, struct Intrusive_List*const faces)
{
	for (struct Intrusive_Link* curr = faces->first; curr; curr = curr->next) {
		struct Face*const face                  = (struct Face*) curr;
		struct Solver_Face_T*const s_face       = (struct Solver_Face_T*) curr;
		struct DG_Solver_Face_T*const dg_s_face = (struct DG_Solver_Face_T*) curr;

		const struct const_Vector_R jdet_fc    = interpret_const_Multiarray_as_Vector_R(s_face->jacobian_det_fc);
		const struct const_Matrix_R normals_fc = interpret_const_Multiarray_as_Matrix_R(s_face->normals_fc);

		struct Matrix_R*const jdet_n_fc = constructor_mm_diag_Matrix_R(1.0,&normals_fc,&jdet_fc,'L',false); // dest.
		if (jdet_n_fc->layout != 'C')
			transpose_Matrix_R(jdet_n_fc,true);

		const struct Test_Case_T*const test_case = (struct Test_Case_T*) sim->test_case_rc->tc;
		const int ind_num_flux = test_case->ind_num_flux[1];
		if (!face->boundary) {
			const bool coll = sim->collocated;
			switch (test_case->solver_method_curr) {
			case 'e': {
				const struct const_Matrix_T*const diff_s_num_s =
					constructor_diff_s_num_s(ind_num_flux,dg_s_face,sim); // destructed

				const struct const_Matrix_T*const*const jdet_n_diff_fc =
					constructor_jdet_n_diff_fc((struct const_Matrix_R*)jdet_n_fc,diff_s_num_s); // destructed
				destructor_const_Matrix_T(diff_s_num_s);

				compute_g_coef_f_i(0,jdet_n_diff_fc,dg_s_face,coll);

				/* Assuming that the numerical solution is central, s_num - s_r == -(s_num - s_l) with
				 * permutation. Noting that the normal is also negated when viewing from the opposite volume,
				 * the "-ve" sign from the normal vector is used to cancel the "-ve" sign above and thus neither
				 * negation is added below. All that is required is then to reorder the jacobian, normal and
				 * numerical solution difference terms as if they are seen from the opposite volume. */
				assert_numerical_solution_is_central(ind_num_flux);
				const struct const_Vector_i*const nc_fc = get_operator__nc_fc_T(0,s_face);
				for (int d = 0; d < DIM; ++d)
					permute_rows_Matrix_T_V((struct Matrix_T*)jdet_n_diff_fc[d],nc_fc);

				compute_g_coef_f_i(1,jdet_n_diff_fc,dg_s_face,coll);

				destructor_jdet_n_diff_fc(jdet_n_diff_fc);
				break;
			} case 'i': {
				/* The procedure used below is inefficient for explicit solver methods. One matrix-matrix
				 * multiplication is removed in the computation of grad_coef_f in the implementation for case
				 * 'e'. */
				compute_d_g_coef_f__d_s_coef_i(0,ind_num_flux,(struct const_Matrix_R*)jdet_n_fc,dg_s_face,coll);

				permute_rows_Matrix_R_V(jdet_n_fc,get_operator__nc_fc_T(0,s_face));
				scale_Matrix_R(jdet_n_fc,-1.0);

				compute_d_g_coef_f__d_s_coef_i(1,ind_num_flux,(struct const_Matrix_R*)jdet_n_fc,dg_s_face,coll);

				compute_g_coef_f_i_using_lin(0,dg_s_face);
				compute_g_coef_f_i_using_lin(1,dg_s_face);
				break;
			} default:
				EXIT_ERROR("Unsupported: %c.",test_case->solver_method_curr);
				break;
			}
		} else {
			compute_g_coef_related_boundary(ind_num_flux,(struct const_Matrix_R*)jdet_n_fc,dg_s_face,sim);
		}
		add_face_grad_coef_f_to_volumes(dg_s_face);
		destructor_Matrix_R(jdet_n_fc);
	}
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for a matrix holding the values of \ref DG_Solver_Volume_T::m_inv multiplied by the appropriate
 *         \ref Solver_Element::tw0_vt_fc operator.
 *  \return See brief. */
static const struct const_Matrix_R* constructor_m_inv_tw0_vt_fc
	(const int side_index,                    ///< The side index of the volume associated with the face.
	 struct DG_Solver_Face_T*const dg_s_face, ///< \ref DG_Solver_Face_T.
	 const bool collocated                    ///< \ref Simulation::collocated.
	);

/** \brief Return a statically allocated array of \ref const_Vector_R\*s with each holding the data in one column of the
 *         input matrix.
 *  \return See brief. */
static struct const_Vector_R* get_jn_fc_V
	(const struct const_Matrix_R*const jdet_n_fc ///< The input matrix.
	);

/** \brief Compute the scaling associated with the \ref Solver_Volume_T::sol_coef contribution to the current face weak
 *         gradient term in the integrated by parts twice form. */
static Real get_sol_scale
	(const int side_index,      ///< The side index of the volume associated with the face.
	 const int sol_index,       ///< The side index of the solution.
	 const int ind_num_flux_2nd ///< The second component of \ref Test_Case_T::ind_num_flux.
	);

static struct Multiarray_Operator get_operator__cv1_vs_vc (const struct DG_Solver_Volume_T*const dg_s_vol)
{
	const struct Volume*const vol               = (struct Volume*) dg_s_vol;
	const struct Solver_Volume_T*const s_vol    = (struct Solver_Volume_T*) dg_s_vol;
	const struct DG_Solver_Element*const dg_s_e = (struct DG_Solver_Element*) vol->element;

	const int p      = s_vol->p_ref,
	          curved = vol->curved;

	return set_MO_from_MO(dg_s_e->cv1_vs_vc[curved],1,(ptrdiff_t[]){0,0,p,p});
}

static const struct const_Matrix_R* constructor_grad_xyz_p
	(const int dir, const struct Multiarray_Operator*const grad_rst, const struct const_Multiarray_R*const metrics)
{
	assert(metrics->layout == 'C');
	const struct const_Matrix_R*const grad_r = grad_rst->data[0]->op_std;
	struct Matrix_R*const grad_xyz = constructor_zero_Matrix_R(grad_r->layout,grad_r->ext_0,grad_r->ext_1); // rtrnd.

	for (int d = 0 ; d < DIM; ++d) {
		const struct const_Vector_R metrics_V =
			{ .ext_0     = metrics->extents[0],
			  .owns_data = false,
			  .data = get_col_const_Multiarray_R(dir*DIM+d,metrics), };
		mm_diag_R('L',1.0,1.0,grad_rst->data[d]->op_std,&metrics_V,grad_xyz,false);
	}

	return (struct const_Matrix_R*) grad_xyz;
}


static void assert_numerical_solution_is_central (const int ind_num_flux_2nd)
{
	switch (ind_num_flux_2nd) {
	case NUM_FLUX_BR2_STABLE: // fallthrough
	case NUM_FLUX_CDG2:
		break; // OK (do nothing)
	default:
		EXIT_ERROR("Central numerical solution **required** by the current implementation. (input: %d).",
		           ind_num_flux_2nd);
		break;
	}
}

static const struct const_Matrix_T* constructor_diff_s_num_s
	(const int ind_num_flux_2nd, const struct DG_Solver_Face_T*const dg_s_face, const struct Simulation*const sim)
{
	const int side_index = 0;
	const struct Solver_Face_T*const s_face = (struct Solver_Face_T*) dg_s_face;

	struct Numerical_Flux_Input_T*const num_flux_i = constructor_Numerical_Flux_Input_T(sim); // destructed
	constructor_Numerical_Flux_Input_data_T(num_flux_i,s_face,sim); // destructed

	const struct Boundary_Value_Input_T*const bv_i = &num_flux_i->bv_l;
	const struct Boundary_Value_T*const bv         = &num_flux_i->bv_r;

	const Real scale[] = { get_sol_scale(side_index,0,ind_num_flux_2nd),
	                       get_sol_scale(side_index,1,ind_num_flux_2nd), };

	const struct const_Multiarray_T*const diff =
		constructor_sum_Multiarrays_const_Multiarray_T(scale[0],bv_i->s,scale[1],bv->s); // destructed/moved

	destructor_Numerical_Flux_Input_data_T(num_flux_i);
	destructor_Numerical_Flux_Input_T(num_flux_i);

	assert(diff->order == 2);
	const struct const_Matrix_T* diff_s_num_s =
		constructor_move_const_Matrix_T_T(diff->layout,diff->extents[0],diff->extents[1],true,diff->data); // rtrnd

	const_cast_b(&diff->owns_data,false);
	destructor_const_Multiarray_T(diff);

	return diff_s_num_s;
}

static const struct const_Matrix_T*const* constructor_jdet_n_diff_fc
	(const struct const_Matrix_R*const jdet_n_fc, const struct const_Matrix_T*const diff_s_num_s)
{
	const struct const_Vector_R*const jn_fc_V = get_jn_fc_V(jdet_n_fc);

	const struct const_Matrix_T** jdet_n_diff_fc = malloc(DIM * sizeof *jdet_n_diff_fc); // free
	for (int d = 0; d < DIM; ++d)
		jdet_n_diff_fc[d] = constructor_mm_diag_const_Matrix_T_R(1.0,diff_s_num_s,&jn_fc_V[d],'L',false); // dest.

	return jdet_n_diff_fc;
}

static void destructor_jdet_n_diff_fc (const struct const_Matrix_T*const* jdet_n_diff_fc)
{
	for (int d = 0; d < DIM; ++d)
		destructor_const_Matrix_T(jdet_n_diff_fc[d]);
	free((void*)jdet_n_diff_fc);
}

static void compute_g_coef_f_i
	(const int side_index, const struct const_Matrix_T*const*const jdet_n_diff_fc,
	 struct DG_Solver_Face_T*const dg_s_face, const bool collocated)
{
	struct Neigh_Info_DG*const ni = &dg_s_face->neigh_info[side_index];

	const struct const_Matrix_R*const m_inv_tw0_vt_fc =
		constructor_m_inv_tw0_vt_fc(side_index,dg_s_face,collocated); // destructed

	for (int d = 0; d < DIM; ++d) {
		struct Multiarray_T grad_coef_f = interpret_Multiarray_as_slice_T(ni->grad_coef_f,2,(ptrdiff_t[]){d});
		struct Matrix_T grad_coef_f_M = interpret_Multiarray_as_Matrix_T(&grad_coef_f);
		mm_RTT('N','N',1.0,0.0,m_inv_tw0_vt_fc,jdet_n_diff_fc[d],&grad_coef_f_M);
	}

	destructor_const_Matrix_R(m_inv_tw0_vt_fc);
}

static void compute_d_g_coef_f__d_s_coef_i
	(const int side_index, const int ind_num_flux_2nd, const struct const_Matrix_R*const jdet_n_fc,
	 struct DG_Solver_Face_T*const dg_s_face, const bool collocated)
{
	const struct Face*const face            = (struct Face*) dg_s_face;
	const struct Solver_Face_T*const s_face = (struct Solver_Face_T*) dg_s_face;

	assert(!face->boundary);

	struct Neigh_Info_DG*const ni = &dg_s_face->neigh_info[side_index];

	const struct const_Matrix_R*const m_inv_tw0_vt_fc =
		constructor_m_inv_tw0_vt_fc(side_index,dg_s_face,collocated); // destructed
	const struct const_Vector_R*const jn_fc_V = get_jn_fc_V(jdet_n_fc);

	for (int sol_index = 0; sol_index < 2; ++sol_index) {
		const struct Operator* cv0_vs_fc_op = get_operator__cv0_vs_fc_T(sol_index,s_face);
		const struct const_Matrix_R* cv0_vs_fc = NULL;
		if (side_index == sol_index) {
			cv0_vs_fc = cv0_vs_fc_op->op_std;
		} else {
			const struct const_Vector_i* nc_fc = get_operator__nc_fc_T(sol_index,s_face);
			cv0_vs_fc = constructor_copy_permute_const_Matrix_R(cv0_vs_fc_op->op_std,nc_fc,'R'); // destructed
		}
		struct Matrix_R*const right = constructor_empty_Matrix_R('R',cv0_vs_fc->ext_0,cv0_vs_fc->ext_1); // dest.

		const Real sol_scale = get_sol_scale(side_index,sol_index,ind_num_flux_2nd);

		for (int d = 0; d < DIM; ++d) {
			mm_diag_R('L',sol_scale,0.0,cv0_vs_fc,&jn_fc_V[d],right,false);

			destructor_const_Matrix_R(ni->d_g_coef_f__d_s_coef[sol_index][d]);
			ni->d_g_coef_f__d_s_coef[sol_index][d] =
				constructor_mm_const_Matrix_R('N','N',1.0,m_inv_tw0_vt_fc,(struct const_Matrix_R*)right,'R'); // keep
		}
		destructor_Matrix_R(right);

		if (side_index != sol_index)
			destructor_const_Matrix_R(cv0_vs_fc);
	}
	destructor_const_Matrix_R(m_inv_tw0_vt_fc);
}

static void compute_g_coef_f_i_using_lin (const int side_index, struct DG_Solver_Face_T*const dg_s_face)
{
	const struct Face*const face = (struct Face*) dg_s_face;
	const struct Solver_Volume_T*const s_vol[] = { (struct Solver_Volume_T*) face->neigh_info[0].volume,
	                                               (struct Solver_Volume_T*) face->neigh_info[1].volume, };
	assert(!face->boundary);

	const struct const_Multiarray_T*const sol_coef[] = { (struct const_Multiarray_T*) s_vol[0]->sol_coef,
	                                                     (struct const_Multiarray_T*) s_vol[1]->sol_coef, };

	struct Neigh_Info_DG*const ni = &dg_s_face->neigh_info[side_index];
	for (int d = 0; d < DIM; ++d) {
		struct Multiarray_T grad_coef_f = interpret_Multiarray_as_slice_T(ni->grad_coef_f,2,(ptrdiff_t[]){d});
		mm_NNC_Multiarray_T(1.0,0.0,ni->d_g_coef_f__d_s_coef[0][d],sol_coef[0],&grad_coef_f);
		mm_NNC_Multiarray_T(1.0,1.0,ni->d_g_coef_f__d_s_coef[1][d],sol_coef[1],&grad_coef_f);
	}
}

static void add_face_grad_coef_f_to_volumes (const struct DG_Solver_Face_T*const dg_s_face)
{
	const struct Face*const face = (struct Face*) dg_s_face;
	const struct Solver_Volume_T*const s_vol[] = { (struct Solver_Volume_T*) face->neigh_info[0].volume,
	                                               (struct Solver_Volume_T*) face->neigh_info[1].volume, };

	const struct Neigh_Info_DG*const ni = dg_s_face->neigh_info;

	const int n_side = ( face->boundary ? 1 : 2 );
	for (int i = 0; i < n_side; ++i)
		add_in_place_Multiarray_T(1.0,s_vol[i]->grad_coef,(struct const_Multiarray_T*)ni[i].grad_coef_f);
}

static void compute_g_coef_related_boundary
	(const int ind_num_flux_2nd, const struct const_Matrix_R*const jdet_n_fc, struct DG_Solver_Face_T*const dg_s_face,
	 const struct Simulation*const sim)
{
	const int side_index = 0,
	          sol_index  = 0;

	const struct Solver_Face_T*const s_face = (struct Solver_Face_T*) dg_s_face;

	const struct Operator*const tw0_vt_fc_op = get_operator__tw0_vt_fc_T(side_index,s_face);
	const struct Operator* cv0_vs_fc_op      = get_operator__cv0_vs_fc_T(sol_index,s_face);

	const struct const_Matrix_R*const tw0_vt_fc = tw0_vt_fc_op->op_std;
	const struct const_Matrix_R* cv0_vs_fc      = cv0_vs_fc_op->op_std;
	struct Neigh_Info_DG*const neigh_info = &dg_s_face->neigh_info[side_index];

	const struct Test_Case_T*const test_case = (struct Test_Case_T*) sim->test_case_rc->tc;
	const int n_vr = test_case->n_var;

	const ptrdiff_t ext_0 = tw0_vt_fc->ext_0,
	                ext_1 = cv0_vs_fc->ext_1,
			    n_fc  = tw0_vt_fc->ext_1;
	assert(n_fc == cv0_vs_fc->ext_0);


	const bool collocated = sim->collocated;
	const struct const_Matrix_R*const m_inv_tw0_vt_fc =
		constructor_m_inv_tw0_vt_fc(side_index,dg_s_face,collocated); // destructed
	const struct const_Vector_R*const jn_fc_V = get_jn_fc_V(jdet_n_fc);

	struct Numerical_Flux_Input_T*const num_flux_i = constructor_Numerical_Flux_Input_T(sim); // destructed
	constructor_Numerical_Flux_Input_data_T(num_flux_i,s_face,sim); // destructed

	const struct Boundary_Value_Input_T*const bv_i = &num_flux_i->bv_l;
	const struct Boundary_Value_T*const bv         = &num_flux_i->bv_r;

	const Real scale[] = { get_sol_scale(side_index,0,ind_num_flux_2nd),
	                       get_sol_scale(side_index,1,ind_num_flux_2nd), };

	// Standard term
	const struct const_Multiarray_T*const diff_s_num_s =
		constructor_sum_Multiarrays_const_Multiarray_T(scale[0],bv_i->s,scale[1],bv->s); // destructed
	const struct const_Matrix_T diff_s_num_s_M = interpret_const_Multiarray_as_Matrix_T(diff_s_num_s);
	for (int d = 0; d < DIM; ++d) {
		const struct const_Matrix_T*const right =
			constructor_mm_diag_const_Matrix_T_R(1.0,&diff_s_num_s_M,&jn_fc_V[d],'L',false); // destructed

		struct Multiarray_T grad_coef_f =
			interpret_Multiarray_as_slice_T(neigh_info->grad_coef_f,2,(ptrdiff_t[]){d});
		struct Matrix_T grad_coef_f_M = interpret_Multiarray_as_Matrix_T(&grad_coef_f);
		mm_RTT('N','N',1.0,0.0,m_inv_tw0_vt_fc,right,&grad_coef_f_M);
		destructor_const_Matrix_T(right);
	}
	destructor_const_Multiarray_T(diff_s_num_s);

	if (test_case->solver_method_curr == 'i') { // Jacobian of standard term.
		for (int d = 0; d < DIM; ++d) {
			destructor_const_Matrix_R(neigh_info->d_g_coef_f__d_s_coef[sol_index][d]);
			neigh_info->d_g_coef_f__d_s_coef[sol_index][d] =
				constructor_empty_const_Matrix_R('R',n_vr*ext_0,n_vr*ext_1); // keep
		}

		const struct const_Multiarray_T*const ds_ds = bv->ds_ds;
		assert(ds_ds->layout == 'C');

		struct Matrix_T*const local_block = constructor_empty_Matrix_T('R',ext_0,ext_1); // destructed
		struct Matrix_T*const right = constructor_empty_Matrix_T('R',cv0_vs_fc->ext_0,cv0_vs_fc->ext_1); // dest.
		struct Vector_T*const diff_ds_ds_V  = constructor_empty_Vector_T(n_fc); // destructed
		struct Vector_T*const ds_ds_jn_fc_V = constructor_empty_Vector_T(n_fc); // destructed
		for (int vr_i = 0; vr_i < n_vr; ++vr_i) {
		for (int vr_b = 0; vr_b < n_vr; ++vr_b) {
			const int ind_ds_ds = vr_b+n_vr*(vr_i);

			const struct const_Vector_T ds_ds_V =
				{ .ext_0 = ds_ds->extents[0], .data = get_col_const_Multiarray_T(ind_ds_ds,ds_ds), };
			set_to_Vector_Vector_T(diff_ds_ds_V,scale[1],&ds_ds_V);
			add_val_to_Vector_T(diff_ds_ds_V,scale[0]);

			for (int d = 0; d < DIM; ++d) {
				dot_mult_Vector_RT(1.0,&jn_fc_V[d],(struct const_Vector_T*)diff_ds_ds_V,ds_ds_jn_fc_V);

				mm_diag_T('L',1.0,0.0,cv0_vs_fc,(struct const_Vector_T*)ds_ds_jn_fc_V,right,false);

				mm_RTT('N','N',1.0,0.0,m_inv_tw0_vt_fc,(struct const_Matrix_T*)right,local_block);

				set_block_Matrix_R((struct Matrix_R*)neigh_info->d_g_coef_f__d_s_coef[sol_index][d],vr_b*ext_0,
				                   vr_i*ext_1,(struct const_Matrix_R*)local_block,0,0,ext_0,ext_1,'i');
			}
		}}
		destructor_Matrix_T(local_block);
		destructor_Matrix_T(right);
		destructor_Vector_T(diff_ds_ds_V);
		destructor_Vector_T(ds_ds_jn_fc_V);
	}

	destructor_Numerical_Flux_Input_data_T(num_flux_i);
	destructor_Numerical_Flux_Input_T(num_flux_i);

	destructor_const_Matrix_R(m_inv_tw0_vt_fc);
}

// Level 2 ********************************************************************************************************** //

static const struct const_Matrix_R* constructor_m_inv_tw0_vt_fc
	(const int side_index, struct DG_Solver_Face_T*const dg_s_face, const bool collocated)
{
	const struct Face*const face                   = (struct Face*) dg_s_face;
	const struct Solver_Face_T*const s_face        = (struct Solver_Face_T*) dg_s_face;
	const struct Solver_Volume_T*const s_vol       = (struct Solver_Volume_T*) face->neigh_info[side_index].volume;
	const struct DG_Solver_Volume_T*const dg_s_vol = (struct DG_Solver_Volume_T*) s_vol;

	const struct Operator*const tw0_vt_fc = get_operator__tw0_vt_fc_T(side_index,s_face);

	const struct const_Matrix_R* m_inv_tw = NULL;
	if (!collocated) {
		m_inv_tw = constructor_mm_const_Matrix_R('N','N',1.0,dg_s_vol->m_inv,tw0_vt_fc->op_std,'R'); // returned
	} else {
		const struct const_Vector_R*const w_vc  = get_operator__w_vc__s_e_T(s_vol);
		const struct const_Vector_R jdet_vc     = interpret_const_Multiarray_as_Vector_R(s_vol->jacobian_det_vc);
		const struct const_Vector_R*const wJ_vc = constructor_dot_mult_const_Vector_R(1.0,w_vc,&jdet_vc,1); // dest.
		m_inv_tw = constructor_mm_diag_const_Matrix_R(1.0,tw0_vt_fc->op_std,wJ_vc,'L',true); // returned
		destructor_const_Vector_R(wJ_vc);
	}
	return m_inv_tw;
}

static struct const_Vector_R* get_jn_fc_V (const struct const_Matrix_R*const jdet_n_fc)
{
	static struct Vector_R jn_fc_V[DIM];

	for (int d = 0; d < DIM; ++d) {
		jn_fc_V[d].ext_0 = jdet_n_fc->ext_0;
		jn_fc_V[d].data  = (Real*)get_col_const_Matrix_R(d,jdet_n_fc);
	}

	return (struct const_Vector_R*) jn_fc_V;
}

static Real get_sol_scale (const int side_index, const int sol_index, const int ind_num_flux_2nd)
{
	switch (ind_num_flux_2nd) {
	case NUM_FLUX_BR2_STABLE: // fallthrough
	case NUM_FLUX_CDG2:
		// Note: Central numerical solution (u_num == 0.5*(u_L + u_R)).
		return 0.5*( side_index == sol_index ? -1.0 : 1.0 );
		break;
	default:
		EXIT_ERROR("Unsupported: %d",ind_num_flux_2nd);
		break;
	}
}
